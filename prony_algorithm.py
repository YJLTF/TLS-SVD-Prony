"""
TLS-SVD-Prony算法实现
用于风电功率7-15天周期分量分析
"""

import numpy as np
from scipy import signal
from scipy.linalg import svd
import pywt


class TLSPronyAnalyzer:
    """改进的TLS-SVD-Prony算法"""

    def __init__(self):
        self.results = None

    def wavelet_denoise(self, data, wavelet='coif5', level=4):
        """小波去噪"""
        coeffs = pywt.wavedec(data, wavelet, level=level)
        threshold = np.std(coeffs[-1]) * np.sqrt(2 * np.log(len(data)))
        denoised_coeffs = [pywt.threshold(c, threshold, mode='soft') for c in coeffs]
        return pywt.waverec(denoised_coeffs, wavelet)[:len(data)]

    def preprocessing(self, data, target_period_min=7, target_period_max=15, fs=1):
        """数据预处理：降采样、去噪、趋势消除"""
        # 去除均值
        data = data - np.mean(data)

        # 趋势消除（二次拟合去除）
        x = np.arange(len(data))
        coeffs = np.polyfit(x, data, 2)
        trend = np.polyval(coeffs, x)
        data = data - trend

        return data

    def tls_prony(self, data, order=None):
        """总体最小二乘Prony算法"""
        n = len(data)

        if order is None:
            order = n // 3

        # 构建Hankel矩阵
        m = n - order
        H = np.zeros((m, order))
        for i in range(m):
            H[i, :] = data[i:i+order]

        # SVD分解
        U, s, Vt = svd(H, full_matrices=False)

        # 确定有效阶数（基于奇异值阈值）
        threshold = s[0] * 1e-6
        effective_order = np.sum(s > threshold)

        if effective_order < 2:
            effective_order = 2

        # TLS降阶
        V = Vt.T
        V_reduced = V[:, :effective_order]

        # 构建Toeplitz矩阵求解特征值
        Z = np.zeros((effective_order, effective_order))
        for i in range(effective_order - 1):
            Z[i, i+1] = 1

        # 最后一行由TLS求解
        try:
            Z[-1, :] = -np.linalg.solve(V_reduced[:effective_order-1, :].T @ V_reduced[:effective_order-1, :],
                                        V_reduced[:effective_order-1, :].T @ V_reduced[effective_order-1, :])
        except:
            Z[-1, :] = 0

        # 求解特征值
        eigenvalues, eigenvectors = np.linalg.eig(Z)

        # 计算频率和阻尼比
        dt = 1.0
        omega = np.log(np.abs(eigenvalues)) / dt
        damping_ratio = np.real(omega) / np.abs(np.log(np.abs(eigenvalues)) + 1e-10)

        # 计算幅值和相位
        amplitudes = np.zeros(len(eigenvalues), dtype=complex)
        for i in range(len(eigenvalues)):
            if i < len(eigenvectors):
                try:
                    amplitudes[i] = np.abs(np.sum(eigenvectors[:, i]))
                except:
                    amplitudes[i] = 0

        # 计算周期
        periods = 2 * np.pi / np.abs(np.imag(eigenvalues))
        periods = np.where(periods > 0, periods, np.inf)

        # 过滤无效结果
        valid_mask = (periods > 0) & (periods < 1e6) & (np.abs(amplitudes) > 0)

        results = []
        for i in range(len(eigenvalues)):
            if valid_mask[i]:
                results.append({
                    'frequency': np.abs(np.imag(eigenvalues[i])) / (2 * np.pi),
                    'period': periods[i],
                    'damping_ratio': damping_ratio[i] if i < len(damping_ratio) else 0,
                    'amplitude': np.abs(amplitudes[i]),
                    'phase': np.angle(eigenvalues[i])
                })

        # 按幅值排序
        results = sorted(results, key=lambda x: x['amplitude'], reverse=True)

        return results

    def calculate_ssnr(self, data, reconstructed):
        """计算信号子空间信噪比"""
        signal_power = np.mean(np.abs(data) ** 2)
        noise_power = np.mean(np.abs(data - reconstructed) ** 2)
        if noise_power > 0:
            return 10 * np.log10(signal_power / noise_power)
        return float('inf')

    def analyze(self, data, params=None):
        """
        完整分析流程

        参数:
            data: 输入数据（一维numpy数组）
            params: 参数字典
                - wavelet: 小波基类型 (默认 'coif5')
                - level: 分解层数 (默认 4)
                - target_period_min: 目标周期最小值（天）(默认 7)
                - target_period_max: 目标周期最大值（天）(默认 15)
                - sampling_interval: 采样间隔（小时）(默认 1)
        """
        if params is None:
            params = {}

        wavelet = params.get('wavelet', 'coif5')
        level = params.get('level', 4)
        target_period_min = params.get('target_period_min', 7)
        target_period_max = params.get('target_period_max', 15)
        sampling_interval = params.get('sampling_interval', 1)

        # 数据预处理
        processed_data = self.preprocessing(data, target_period_min, target_period_max)

        # 小波去噪
        denoised_data = self.wavelet_denoise(processed_data, wavelet, level)

        # TLS-Prony分析
        prony_results = self.tls_prony(denoised_data)

        # 筛选7-15天周期分量
        target_periods = []
        for r in prony_results:
            period_days = r['period'] * sampling_interval / 24
            r['period_days'] = period_days
            if target_period_min <= period_days <= target_period_max:
                target_periods.append(r)

        # 存储结果
        self.results = {
            'all_components': prony_results,
            'target_components': target_periods,
            'denoised_data': denoised_data,
            'original_data': data,
            'params': params
        }

        return self.results

    def get_periodogram(self, data):
        """计算功率谱密度"""
        frequencies, psd = signal.periodogram(data, fs=1.0)
        periods = 1 / frequencies[1:]
        psd_db = 10 * np.log10(psd[1:] + 1e-10)
        return periods, psd_db, frequencies

    def reconstruct_signal(self, components, n_points):
        """重构信号"""
        t = np.arange(n_points)
        reconstructed = np.zeros(n_points, dtype=complex)

        for comp in components:
            omega = 2 * np.pi * comp['frequency']
            damping = comp['damping_ratio']
            amplitude = comp['amplitude']
            phase = comp['phase']

            reconstructed += amplitude * np.exp((damping + 1j * omega) * t + 1j * phase)

        return np.real(reconstructed)


def generate_sample_data(n_days=30, sampling_interval=1):
    """生成示例数据（包含7-15天周期分量）"""
    n_points = int(n_days * 24 / sampling_interval)
    t = np.arange(n_points) * sampling_interval / 24  # 转换为天

    # 信号组成
    # 10天周期分量
    signal_10d = 2.0 * np.sin(2 * np.pi * t / 10 + 0.5)
    # 12天周期分量
    signal_12d = 1.5 * np.sin(2 * np.pi * t / 12 + 1.2)
    # 随机噪声
    noise = np.random.randn(n_points) * 0.3

    # 组合信号
    data = signal_10d + signal_12d + noise + 5  # 加上直流分量

    return data


if __name__ == '__main__':
    # 测试代码
    analyzer = TLSPronyAnalyzer()

    # 生成示例数据
    data = generate_sample_data(n_days=30)

    # 分析
    results = analyzer.analyze(data)

    print("检测到的周期分量:")
    for comp in results['target_components']:
        print(f"  周期: {comp['period_days']:.2f} 天, "
              f"幅值: {comp['amplitude']:.4f}, "
              f"阻尼比: {comp['damping_ratio']:.4f}")

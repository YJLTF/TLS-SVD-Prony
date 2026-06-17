"""
风电功率周期分析Web应用
基于改进TLS-SVD-Prony算法
"""

import io
import base64
from datetime import datetime
from flask import Flask, render_template, request, jsonify, send_file
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prony_algorithm import TLSPronyAnalyzer, generate_sample_data

app = Flask(__name__)
app.secret_key = 'tls-svd-prony-secret-key-2024'


@app.route('/')
def index():
    """主页"""
    return render_template('index.html')


@app.route('/api/analyze', methods=['POST'])
def analyze():
    """分析数据"""
    try:
        # 获取上传的文件
        file = request.files.get('datafile')
        params = request.form.to_dict()

        # 解析参数
        analysis_params = {
            'wavelet': params.get('wavelet', 'coif5'),
            'level': int(params.get('level', 4)),
            'target_period_min': float(params.get('target_period_min', 7)),
            'target_period_max': float(params.get('target_period_max', 15)),
            'sampling_interval': float(params.get('sampling_interval', 1))
        }

        # 读取数据
        if file:
            content = file.read().decode('utf-8')
            data = np.loadtxt(io.StringIO(content))
        else:
            # 生成示例数据
            n_days = int(params.get('n_days', 30))
            data = generate_sample_data(n_days=n_days, sampling_interval=analysis_params['sampling_interval'])

        # 分析
        analyzer = TLSPronyAnalyzer()
        results = analyzer.analyze(data, analysis_params)

        # 生成图表
        plots = generate_plots(results, analysis_params)

        # 准备返回结果
        response = {
            'success': True,
            'target_components': results['target_components'][:10],  # 最多10个分量
            'all_components_count': len(results['all_components']),
            'target_count': len(results['target_components']),
            'params': analysis_params,
            'plots': plots
        }

        return jsonify(response)

    except Exception as e:
        import traceback
        return jsonify({
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc()
        })


def generate_plots(results, params):
    """生成图表并返回base64编码"""
    plots = {}

    # 1. 时域波形对比图
    fig1, axes = plt.subplots(2, 1, figsize=(12, 8))

    # 原始数据
    axes[0].plot(results['original_data'][:500], 'b-', alpha=0.7, linewidth=0.8)
    axes[0].set_title('原始信号（局部）', fontsize=14)
    axes[0].set_xlabel('采样点')
    axes[0].set_ylabel('幅值')
    axes[0].grid(True, alpha=0.3)

    # 去噪后数据
    axes[1].plot(results['denoised_data'][:500], 'r-', alpha=0.7, linewidth=0.8)
    axes[1].set_title('小波去噪后信号（局部）', fontsize=14)
    axes[1].set_xlabel('采样点')
    axes[1].set_ylabel('幅值')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    buf1 = io.BytesIO()
    plt.savefig(buf1, format='png', dpi=100, bbox_inches='tight')
    buf1.seek(0)
    plots['time_domain'] = base64.b64encode(buf1.read()).decode('utf-8')
    plt.close(fig1)

    # 2. 周期分量分布图
    fig2, ax = plt.subplots(figsize=(12, 6))

    target = results['target_components']
    if target:
        periods = [c['period_days'] for c in target]
        amplitudes = [c['amplitude'] for c in target]

        bars = ax.bar(range(len(periods)), amplitudes, color='steelblue', alpha=0.8)
        ax.set_xticks(range(len(periods)))
        ax.set_xticklabels([f'{p:.1f}天' for p in periods], rotation=45)
        ax.set_title('检测到的7-15天周期分量', fontsize=14)
        ax.set_xlabel('周期（天）')
        ax.set_ylabel('幅值')
        ax.grid(True, alpha=0.3, axis='y')

        # 在柱状图上标注数值
        for bar, amp in zip(bars, amplitudes):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                   f'{amp:.2f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    buf2 = io.BytesIO()
    plt.savefig(buf2, format='png', dpi=100, bbox_inches='tight')
    buf2.seek(0)
    plots['period_components'] = base64.b64encode(buf2.read()).decode('utf-8')
    plt.close(fig2)

    # 3. 功率谱密度图
    fig3, ax = plt.subplots(figsize=(12, 6))

    analyzer = TLSPronyAnalyzer()
    periods, psd_db, freqs = analyzer.get_periodogram(results['denoised_data'])

    # 只显示周期大于1天的部分
    valid = periods >= 1
    ax.plot(periods[valid], psd_db[valid], 'g-', linewidth=1.2)
    ax.axvline(x=7, color='r', linestyle='--', alpha=0.7, label='7天')
    ax.axvline(x=15, color='r', linestyle='--', alpha=0.7, label='15天')
    ax.set_title('功率谱密度（周期域）', fontsize=14)
    ax.set_xlabel('周期（天）')
    ax.set_ylabel('功率谱密度 (dB)')
    ax.set_xlim([0, 50])
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    buf3 = io.BytesIO()
    plt.savefig(buf3, format='png', dpi=100, bbox_inches='tight')
    buf3.seek(0)
    plots['power_spectrum'] = base64.b64encode(buf3.read()).decode('utf-8')
    plt.close(fig3)

    # 4. 重构信号对比图
    fig4, ax = plt.subplots(figsize=(12, 5))

    # 取前200点进行重构对比
    n_show = min(200, len(results['original_data']))
    t = np.arange(n_show)

    # 原始信号
    ax.plot(t, results['original_data'][:n_show], 'b-', alpha=0.5,
            linewidth=1, label='原始信号')

    # 重构信号（只使用目标周期分量）
    if target:
        reconstructed = np.zeros(n_show)
        sampling_interval = params.get('sampling_interval', 1)
        for comp in target:
            omega = 2 * np.pi * comp['frequency']
            damping = comp['damping_ratio']
            amplitude = comp['amplitude']
            phase = comp['phase']
            reconstructed += amplitude * np.exp((damping + 1j * omega) * t * sampling_interval / 24 + 1j * phase)
        ax.plot(t, np.real(reconstructed), 'r-', linewidth=2, label='重构信号(7-15天分量)')

    ax.set_title('原始信号与重构信号对比', fontsize=14)
    ax.set_xlabel('采样点')
    ax.set_ylabel('幅值')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    buf4 = io.BytesIO()
    plt.savefig(buf4, format='png', dpi=100, bbox_inches='tight')
    buf4.seek(0)
    plots['reconstruction'] = base64.b64encode(buf4.read()).decode('utf-8')
    plt.close(fig4)

    return plots


@app.route('/api/download_sample', methods=['GET'])
def download_sample():
    """下载示例数据"""
    data = generate_sample_data(n_days=30)
    output = io.StringIO()
    for value in data:
        output.write(f'{value}\n')
    output.seek(0)

    return send_file(
        io.BytesIO(output.getvalue().encode('utf-8')),
        mimetype='text/csv',
        as_attachment=True,
        download_name='sample_wind_power_data.csv'
    )


if __name__ == '__main__':
    print("启动风电功率周期分析系统...")
    print("访问 http://127.0.0.1:5000 查看应用")
    app.run(debug=True, host='0.0.0.0', port=5000)

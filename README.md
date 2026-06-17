# TLS-SVD-Prony

基于改进TLS-SVD-Prony算法的风电功率7-15天周期分量分析方法

## 项目简介

本项目提供一个基于Web的交互式分析工具，用于从风电功率时间序列数据中提取并验证7-15天周期的分量。采用改进的TLS-SVD-Prony算法，结合小波变换与奇异值分解（SVD），有效提升抗噪能力与参数估计稳定性。

## 核心算法

- **TLS-SVD-Prony** - 总体最小二乘结合奇异值分解的Prony算法
- **小波去噪** - 支持多种小波基（Coif5、Daubechies、Symlet等）
- **SVD定阶** - 基于奇异值阈值自动确定有效阶数
- **周期分量提取** - 精准识别7-15天周期范围

## 目录结构

```
.
├── [app.py](file:///workspace/app.py)              # Flask Web应用
├── [prony_algorithm.py](file:///workspace/prony_algorithm.py)   # TLS-SVD-Prony算法核心实现
├── [requirements.txt](file:///workspace/requirements.txt)      # Python依赖
└── [templates/](file:///workspace/templates/)
    └── [index.html](file:///workspace/templates/index.html)   # 前端UI界面
```

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

### 2. 启动Web服务

```bash
python app.py
```

### 3. 访问应用

在浏览器中打开: http://127.0.0.1:5000

## 使用说明

### 数据上传

- **支持格式**: CSV、TXT、DAT
- **数据要求**: 每行一个数值，代表一个采样点的功率值
- **采样间隔**: 可自定义（默认1小时）

### 参数调整

| 参数 | 说明 | 默认值 |
|------|------|--------|
| 小波基类型 | 选择小波变换的基函数 | Coif5 |
| 分解层数 | 小波分解层数（1-10） | 4 |
| 目标周期最小值 | 期望检测的周期下限（天） | 7 |
| 目标周期最大值 | 期望检测的周期上限（天） | 15 |
| 采样间隔 | 数据采样间隔（小时） | 1 |

### 可视化结果

1. **时域信号对比** - 原始信号与去噪后信号对比
2. **周期分量分布** - 检测到的7-15天周期分量幅值分布
3. **功率谱密度** - 信号在周期域上的频谱分析
4. **信号重构对比** - 使用目标周期分量重构信号与原始信号对比
5. **分量详情表格** - 每个检测到的周期分量的详细参数（周期、频率、幅值、阻尼比、相位）

## 技术栈

- **后端**: Flask 3.0
- **数值计算**: NumPy, SciPy
- **信号处理**: PyWavelets
- **可视化**: Matplotlib
- **前端**: 原生 HTML / CSS / JavaScript

## 依赖列表

```
Flask>=3.0.0
numpy>=1.24.0
scipy>=1.11.0
matplotlib>=3.8.0
pywavelets>=1.9.0
pandas>=2.0.0
Werkzeug>=3.0.0
```

## 使用示例

```python
from prony_algorithm import TLSPronyAnalyzer, generate_sample_data

# 创建分析器实例
analyzer = TLSPronyAnalyzer()

# 生成示例数据（或读取真实数据）
data = generate_sample_data(n_days=30)

# 设置参数并分析
params = {
    'wavelet': 'coif5',
    'level': 4,
    'target_period_min': 7,
    'target_period_max': 15,
    'sampling_interval': 1
}

results = analyzer.analyze(data, params)

# 查看结果
print(f"检测到 {len(results['target_components'])} 个目标周期分量")
for comp in results['target_components']:
    print(f"  周期: {comp['period_days']:.2f} 天, 幅值: {comp['amplitude']:.4f}")
```

## License

MIT License

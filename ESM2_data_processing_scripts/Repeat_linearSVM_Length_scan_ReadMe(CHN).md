# Repeat_linearSVM_Length_scan.py - 植物P450蛋白N端长度扫描分析（多次重复线性SVM版本）

## 脚本功能概述

本脚本用于系统分析N端长度对植物P450蛋白分类性能的影响。采用重复线性SVM方法（默认6次重复），系统扫描40-100个氨基酸的N端序列，验证三个进化假说：快速分化期（40-49AA）、过渡期（50-84AA）和去噪音期（85-100AA）。脚本结合ESM2嵌入向量与稳健的机器学习方法，确保统计结果的可靠性。

## 主要功能模块

1. **重复分析框架**：使用不同随机种子执行6次独立的线性SVM分析，确保结果稳健性
2. **多长度特征提取**：提取不同长度（40-100AA）N端序列的ESM2嵌入向量
3. **假说验证**：专门测试P450 N端进化的三个关键时期
4. **统计整合**：整合6次独立重复结果，计算平均性能、显著比例和效应量
5. **综合可视化**：生成整合图展示性能曲线、显著热图和拐点分析
6. **数据导出**：以CSV和JSON格式导出详细结果，便于进一步分析

## 系统要求与依赖包

### Python版本要求
- Python 3.8 或更高版本

### 必需Python包
```bash
pip install torch transformers numpy pandas matplotlib seaborn scikit-learn tqdm scipy
```

### 高级用户可选安装
```bash
# 可选：GPU加速（需要CUDA兼容的GPU）
# 确保安装合适的CUDA版本
pip install torch --index-url https://download.pytorch.org/whl/cu118  # CUDA 11.8示例
```

## 文件准备

### 必需数据文件
将以下文件放置在脚本同一目录中：

1. **FASTA文件**：包含P450蛋白序列
   - 默认文件名：`P450_unique_pep_final.fasta`
   - 格式：标准FASTA格式，序列头包含基因标识符

2. **聚类标签文件**：包含功能聚类信息
   - 默认文件名：`P450_unique_pep_final_esm_clusters_genes.csv`
   - 必需列：必须包含基因标识符和聚类标签
   - 脚本自动检测分隔符（逗号、制表符、分号、竖线）

### 文件结构示例
```
项目目录/
├── repeat_linearSVM_Length_scan.py      # 主脚本
├── P450_unique_pep_final.fasta          # FASTA序列文件
├── P450_unique_pep_final_esm_clusters_genes.csv  # 聚类标签文件
└── p450_length_scan_repeated_svm/       # 输出目录（运行时创建）
```

## 配置参数

### 命令行参数
脚本支持以下命令行参数：

```python
# 使用默认参数的基本用法
python repeat_linearSVM_Length_scan.py

# 自定义分析
python repeat_linearSVM_Length_scan.py \
    --fasta "my_sequences.fasta" \
    --cluster "my_clusters.csv" \
    --output "custom_output" \
    --min_length 30 \
    --max_length 120 \
    --n_repeats 10 \
    --device "cuda" \
    --seed_start 123
```

### 参数说明
```python
# 默认值及说明：
--fasta: "P450_unique_pep_final.fasta"           # FASTA文件路径
--cluster: "P450_unique_pep_final_esm_clusters_genes.csv"  # 聚类文件路径
--output: "p450_length_scan_repeated_svm"        # 输出目录名称
--min_length: 40                                 # 分析的最小N端长度
--max_length: 100                                # 分析的最大N端长度
--n_repeats: 6                                   # 独立重复次数
--device: "cpu"                                  # 计算设备（"cpu"或"cuda"）
--seed_start: 42                                 # 起始随机种子（每次递增100）
```

## 运行脚本

### 基本执行
```bash
# 1. 确保数据文件在同一个目录
# 2. 使用默认参数运行脚本
python repeat_linearSVM_Length_scan.py
```

### 高级执行
```bash
# 使用自定义参数运行
python repeat_linearSVM_Length_scan.py \
    --min_length 30 \
    --max_length 150 \
    --n_repeats 10 \
    --device "cuda"  # 如有GPU可用
```

### 进度监控
脚本提供详细的进度信息：
- 模型加载状态
- 数据加载统计
- 单个重复进度（1-6）
- 统计整合步骤
- 可视化生成

## 输出结果

### 目录结构
```
p450_length_scan_repeated_svm/                # 主输出目录
├── repeat_1_seed_42/                         # 单次重复结果
│   ├── n_terminal_length_scan.csv
│   └── full_length_performance.csv
├── repeat_2_seed_142/                        # 第二次重复
│   ├── n_terminal_length_scan.csv
│   └── full_length_performance.csv
├── ...                                       # 更多重复（3-6）
├── integrated_results/                       # 整合结果
│   ├── integrated_n_terminal_length_scan.csv
│   ├── integrated_full_length_performance.csv
│   ├── detailed_all_repeats_accuracy.csv
│   ├── all_results_summary.json
│   └── plots/                                # 可视化目录
│       ├── integrated_length_scan_analysis.png/pdf
│       ├── hypothesis_testing_turning_point.png/pdf
│       └── individual_repeats_comparison.png/pdf
└── console_output.txt                        # 运行日志（如重定向）
```

### 关键输出文件

1. **整合分析文件** (`integrated_results/`)
   - `integrated_n_terminal_length_scan.csv`：各长度的平均值和变异度
   - `detailed_all_repeats_accuracy.csv`：每次重复各长度的原始准确率
   - `all_results_summary.json`：综合摘要，包括拐点分析

2. **可视化文件** (`plots/`)
   - `integrated_length_scan_analysis.png/pdf`：四面板图，展示：
     - 平均准确率曲线及标准差
     - 显著比例热图
     - 性能差距百分比
     - 准确率分布箱线图
   - `hypothesis_testing_turning_point.png/pdf`：主要假说验证图
   - `individual_repeats_comparison.png/pdf`：6次重复的对比图

3. **单次重复文件** (`repeat_*/`)
   - `n_terminal_length_scan.csv`：单次重复的完整结果
   - `full_length_performance.csv`：全长序列性能

## 结果解读指南

### 关键统计指标

1. **拐点识别**
   - 脚本识别第一个与全长序列性能无显著差异的N端长度
   - 标准：显著比例<50%，准确率≥全长性能的90%

2. **显著比例**
   - 显示与全长有显著差异（p<0.05）的重复比例
   - 高比例（>0.7）：有显著差异的强证据
   - 低比例（<0.3）：无显著差异的强证据

3. **性能差距**
   - 与全长准确率的百分比差异
   - 负值：N端性能优于全长
   - 正值：N端性能劣于全长

### 假说验证区域

1. **快速分化期（40-49AA）**
   - 预期：高显著比例，大的正性能差距
   - 解释：短N端片段不足以进行准确分类

2. **过渡期（50-84AA）**
   - 预期：显著比例降低，性能差距减小
   - 解释：功能信息逐渐积累

3. **去噪音期（85-100AA）**
   - 预期：低显著比例，接近零或负性能差距
   - 解释：额外序列去除噪音但不增加信息

### 生物学解释示例

```python
# 强N端信号：
拐点 ≤ 55AA，准确率比例 ≥ 0.95
# 解释：功能信息高度集中在N端区域

# 分散信号：
拐点 ≥ 80AA，准确率比例 < 0.90
# 解释：功能信息分布在序列多个区域
```

## 方法学细节

### ESM2模型配置
- 使用 `facebook/esm2_t12_35M_UR50D`（3500万参数）
- 备选模型：`facebook/esm2_t6_8M_UR50D`（如主模型不可用）
- 嵌入提取：所有残基嵌入向量的平均池化

### SVM配置
- 使用LinearSVC，自动参数调优（C值：0.01、0.1、1.0、10.0）
- 类别权重：对不平衡数据集使用平衡权重
- PCA降维（保留90%方差）
- 5折分层交叉验证

### 稳健性特性
- 多配置SVM及后备机制
- 自动处理稀有类别（<5个样本→合并）
- 使用随机SVD的PCA提高数值稳定性
- 内存高效的批处理

## 故障排除

### 常见问题与解决方案

1. **内存问题**
```python
# 解决方案1：使用较小的ESM2模型
# 修改脚本第34行：
model_name="facebook/esm2_t6_8M_UR50D"

# 解决方案2：减小批处理大小
# 修改脚本第271行：
batch_size = 1  # 默认值为2
```

2. **运行缓慢**
```python
# 使用GPU加速
python repeat_linearSVM_Length_scan.py --device "cuda"

# 减少重复次数
python repeat_linearSVM_Length_scan.py --n_repeats 3
```

3. **数据加载错误**
```bash
# 检查文件格式
# FASTA：确保标准格式，>开头为序列头
# CSV：确保基因标识符与FASTA头匹配
```

4. **可视化错误**
```bash
# 安装必要包
pip install seaborn
# 或在代码中禁用可视化（注释第1199-1369行）
```

### 错误信息与解决方案

1. **"模型加载失败"**
   - 检查模型下载的网络连接
   - 尝试较小模型：`facebook/esm2_t6_8M_UR50D`

2. **"未找到共同基因"**
   - 验证FASTA和CSV文件的基因标识符是否匹配
   - 检查CSV文件分隔符检测

3. **"CUDA内存不足"**
   - 使用CPU：`--device "cpu"`
   - 在代码中减小批处理大小（第271行）

## 高级定制

### 修改分析参数
```python
# 在main()函数（第1296-1331行）中修改：
analyzer = RepeatedLinearSVMLengthScanAnalyzer(
    model_name="facebook/esm2_t33_650M_UR50D",  # 更大模型
    device="cuda"  # 强制使用GPU
)

# 在robust_linear_svm_cv方法（第342-347行）中调整：
svm_configs = [
    {'C': 0.001, 'max_iter': 10000},  # 更保守
    {'C': 0.01, 'max_iter': 10000},
    {'C': 0.1, 'max_iter': 10000},
    {'C': 1.0, 'max_iter': 20000},    # 更多迭代
    {'C': 10.0, 'max_iter': 20000},
]
```

### 添加新假说区域
```python
# 在analyze_hypothesis_support方法（第1109-1115行）中添加：
hypothesis_regions = {
    'rapid_divergence': (40, 49),
    'transition': (50, 84),
    'noise_removal': (85, 100),
    'new_hypothesis': (101, 120)  # 自定义区域
}
```

### 扩展长度范围
```python
# 命令行选项：
python repeat_linearSVM_Length_scan.py --min_length 20 --max_length 200

# 或修改参数解析器的默认值（第1273-1276行）
parser.add_argument('--min_length', type=int, default=20)
parser.add_argument('--max_length', type=int, default=200)
```

## 引用与致谢

如使用本脚本，请引用：

- **ESM2模型**：Lin et al. (2022) *Nature Methods*
- **线性SVM**：Scikit-learn库 (Pedregosa et al., 2011)
- **统计方法**：SciPy库 (Virtanen et al., 2020)

**方法学参考文献**：  
*植物P450进化中N端功能信号的稳健检测：基于重复线性SVM的分析方法*

## 技术支持与联系

技术问题：
1. 检查控制台错误信息
2. 验证数据文件格式和位置
3. 确保所有依赖包已安装
4. 减小分析规模进行调试

关于方法学或结果解读的科学问题，请参考随附论文或联系通讯作者。

---

**脚本版本**：repeat_linearSVM_Length_scan.py v1.0  
**最后更新**：2025年  
**设计目的**：使用重复线性SVM分析对P450 N端进化进行稳健假说验证
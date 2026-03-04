# Classifier_Selection.py - 植物P450蛋白区域分类性能比较分析脚本

## 脚本功能概述

该脚本用于系统比较植物P450蛋白不同序列区域（N端、C端、中间区域、随机片段、全长序列）在功能分类任务中的性能差异。通过结合ESM2蛋白语言模型、多种机器学习分类器和降维可视化方法，本脚本旨在验证N端序列信号的特异性，区分序列顺序信息与氨基酸组成信息的贡献，为理解P450蛋白功能进化的结构基础提供定量依据。

## 主要功能模块

1. **多区域特征提取**：提取N端、C端、中间区域、随机片段和全长序列的ESM2嵌入向量
2. **对照实验设计**：
   - Shuffled N端：打乱N端氨基酸顺序，控制氨基酸组成
   - 随机片段：从序列随机位置提取片段
   - 多次随机片段试验：评估随机位置的稳健性
3. **分类器比较**：比较LinearSVM、LogisticRegression和RandomForest的分类性能
4. **多维可视化**：提供PCA、t-SNE、UMAP三种降维方法的结果
5. **统计验证**：包括Bootstrap置信区间、配对t检验和效应量分析
6. **数据导出**：导出所有可视化所需原始数据，便于R语言二次分析

## 系统要求与依赖包

### Python版本要求
- Python 3.8 或更高版本

### 必需Python包
```bash
pip install torch transformers numpy pandas matplotlib seaborn scikit-learn tqdm scipy
```

### 可选Python包（用于高级降维）
```bash
# 用于t-SNE（已在scikit-learn中）
# 用于UMAP（如需要请安装）
pip install umap-learn
```

## 文件准备

### 必需数据文件
将以下文件放置在脚本同一目录中：

1. **FASTA文件**：包含P450蛋白序列
   - 默认文件名：`P450_unique_pep_final.fasta`
   - 格式：标准FASTA格式，序列头包含基因标识符

2. **聚类标签文件**：包含基因到功能类别/亚家族的映射
   - 默认文件名：`P450_unique_pep_final_esm_clusters_genes.csv`
   - 必需列：`Gene_List`（基因标识符列表，分号分隔）和`Cluster`（类别标签）

### 文件结构示例
```
项目目录/
├── Classifier_Selection.py          # 主脚本
├── P450_unique_pep_final.fasta      # FASTA序列文件
├── P450_unique_pep_final_esm_clusters_genes.csv  # 聚类标签文件
└── enhanced_p450_region_analysis/   # 运行后生成的输出目录
```

## 配置参数

### 关键参数说明（位于`main()`函数中）
```python
# 文件路径配置
FASTA_FILE = "P450_unique_pep_final.fasta"
CLUSTER_FILE = "P450_unique_pep_final_esm_clusters_genes.csv"
OUTPUT_DIR = "enhanced_p450_region_analysis"

# 分析器初始化参数
analyzer = EnhancedP450Analyzer(
    model_name="facebook/esm2_t33_650M_UR50D",  # ESM2模型版本
    n_aa=100,      # N端氨基酸分析长度
    c_aa=100,      # C端氨基酸分析长度
    middle_aa=100, # 中间区域分析长度
    random_aa=100, # 随机片段长度
    device="cpu",  # 使用CPU，如有GPU可改为"cuda"
    use_tsne=True, # 启用t-SNE降维
    use_umap=True, # 启用UMAP降维
    num_random_trials=20  # 随机片段重复试验次数
)
```

## 运行脚本

1. **确保数据文件就位**：将FASTA文件和聚类文件放在脚本同目录

2. **运行分析**：
```bash
python Classifier_Selection.py
```

3. **监控输出**：脚本执行过程中会显示进度信息和中间结果

## 输出结果

### 目录结构
```
enhanced_p450_region_analysis/
├── plots/                                      # 可视化图表目录
│   ├── classifier_performance_comparison_bar.png/pdf
│   ├── region_performance_comparison_bar.png/pdf
│   ├── performance_heatmap.png/pdf
│   ├── n_terminal_vs_shuffled_comparison.png/pdf
│   ├── region_dimension_reduction.png/pdf
│   └── random_fragment_distributions.png/pdf
├── *.csv                                      # 原始数据文件
│   ├── N_terminal_projection_data.csv
│   ├── Shuffled_N_terminal_projection_data.csv
│   ├── C_terminal_projection_data.csv
│   ├── Middle_Region_projection_data.csv
│   ├── Random_Fragment_projection_data.csv
│   ├── Full_length_projection_data.csv
│   ├── classifier_performance.csv
│   ├── detailed_accuracy_data.csv
│   ├── bootstrap_distribution.csv
│   ├── class_distribution.csv
│   ├── region_performance_comparison.csv
│   └── summary_statistics.csv
└── 控制台输出统计结果和关键结论
```

### 主要输出文件说明

1. **投影数据文件** (`*_projection_data.csv`)
   - 包含每个区域的PCA/t-SNE/UMAP二维坐标
   - 用于R语言绘制高级定制可视化图表

2. **分类器性能文件** (`classifier_performance.csv`)
   - 每个区域、每个分类器的交叉验证准确率统计

3. **详细准确率数据** (`detailed_accuracy_data.csv`)
   - 每个交叉验证折叠的准确率，适用于方差分析

4. **区域性能比较** (`region_performance_comparison.csv`)
   - 各区域在相同交叉验证划分下的直接比较

5. **统计摘要文件** (`summary_statistics.csv`)
   - 用于制作柱状图的汇总统计数据

## 结果解读指南

### 关键统计指标

1. **N端与Shuffled N端比较**
   - N端 > Shuffled N端：表明存在序列顺序信息
   - N端 ≈ Shuffled N端：仅氨基酸组成信息起作用

2. **N端与其他区域比较**
   - N端性能最优：支持N端信号特异性假说
   - 多区域性能相近：功能信息分布广泛

3. **多次随机片段试验**
   - 随机片段性能稳定：表明分类任务鲁棒性
   - 随机片段性能波动大：结果可能受随机位置影响

### 生物学意义判断标准

1. **强特异性信号**：
   ```
   N端 > Shuffled N端
   N端 > C端
   N端 > 中间区域
   N端 > 随机片段
   ```

2. **序列顺序依赖**：
   ```
   N端 > Shuffled N端
   Cohen's d > 0.5 (中等以上效应量)
   p-value < 0.05 (统计显著)
   ```

## 故障排除

### 常见问题与解决方案

1. **内存不足**
   - 降低批处理大小：修改`batch_size`参数
   - 使用较小ESM2模型：如`facebook/esm2_t12_35M_UR50D`
   - 启用PCA降维：减少特征维度

2. **数据加载失败**
   - 检查FASTA文件格式：确保为标准FASTA格式
   - 检查CSV分隔符：脚本自动检测，也可手动指定
   - 验证基因标识符一致性：确保FASTA和CSV文件中的基因名匹配

3. **模型下载缓慢或失败**
   - 手动下载模型：从Hugging Face下载后指定本地路径
   - 使用代理：设置网络代理
   - 选择较小模型：如`facebook/esm2_t6_8M_UR50D`

4. **UMAP/t-SNE不可用**
   - 安装缺失包：`pip install umap-learn`
   - 禁用对应功能：设置`use_umap=False`或`use_tsne=False`

## 高级定制

### 修改分析区域长度
```python
# 在main()函数中修改
analyzer = EnhancedP450Analyzer(
    n_aa=50,      # 分析N端前50个氨基酸
    c_aa=50,      # 分析C端后50个氨基酸
    # ... 其他参数
)
```

### 增加分类器
```python
# 在compare_feature_regions方法中添加分类器
classifiers = {
    'LinearSVM': LinearSVC(...),
    'LogisticRegression': LogisticRegression(...),
    'RandomForest': RandomForest(...),
    '新增分类器': 分类器实例
}
```

### 调整统计阈值
```python
# 修改统计显著性阈值
alpha = 0.01  # 更严格的显著性水平
# 在paired t-test部分修改
```

## 引用与致谢

如使用本脚本，请引用相关方法：
- ESM2蛋白语言模型：Lin et al. (2022) Nature Methods
- 机器学习框架：scikit-learn
- 降维方法：PCA (scikit-learn), t-SNE (van der Maaten & Hinton, 2008), UMAP (McInnes et al., 2018)

## 技术支持与反馈

如有问题或建议，请：
1. 检查控制台错误信息
2. 查看输出目录中的日志文件
3. 确保所有依赖包版本兼容
4. 提供可重现的示例数据用于调试

---
*脚本设计者：P450进化分析项目组*
*最后更新：2025年*
*版本：EnhancedP450Analyzer v1.0*
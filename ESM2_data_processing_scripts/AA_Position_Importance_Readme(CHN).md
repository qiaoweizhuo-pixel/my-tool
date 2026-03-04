# AA_Position_Importance.py - 植物P450蛋白氨基酸位置重要性分析（6次重复增强版）

## 脚本功能概述

本脚本用于系统分析植物P450蛋白N端序列中各氨基酸位置对功能分类的重要性。采用6次重复的线性SVM方法，结合ESM2嵌入向量，识别在多次独立实验中一致重要的氨基酸位置。脚本还集成了与正选择位点的比较分析，验证功能重要位置与适应性进化位置的关联性。该脚本修复了JSON序列化问题，确保所有结果可被正确保存和读取。

## 主要功能模块

1. **6次重复分析框架**：使用不同随机种子执行6次独立实验，确保位置重要性评估的稳健性
2. **位置特异性特征提取**：提取每个氨基酸位置的独立ESM2嵌入向量
3. **SVM权重反投影**：将PCA空间中的SVM权重反投影到原始特征空间，计算位置重要性
4. **稳定性评估**：计算位置重要性的变异系数、权重相关性、Kendall's W一致性系数
5. **正选择位点整合**：与用户提供的正选择位点进行超几何检验和Jaccard相似度分析
6. **综合可视化**：生成位置重要性条形图、热图、累积分布图、稳定性分析图
7. **详细报告生成**：自动生成包含统计分析、生物学解释和实验建议的综合报告

## 系统要求与依赖包

### Python版本要求

- Python 3.8 或更高版本

### 必需Python包

```bash
pip install torch transformers numpy pandas matplotlib seaborn scikit-learn tqdm scipy
```

### 可选的GPU支持

```bash
# 如需GPU加速，安装对应CUDA版本的PyTorch
# 示例：CUDA 11.8
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

## 文件准备

### 必需数据文件

将以下文件放置在脚本同一目录中：

1. **FASTA文件**：包含P450蛋白序列
   - 默认文件名：`P450_unique_pep_final.fasta`
   - 格式：标准FASTA格式，序列头包含基因标识符

2. **聚类标签文件**：包含功能聚类信息
   - 默认文件名：`P450_unique_pep_final_esm_clusters_genes.csv`
   - 必需列：`Gene_List`（基因列表，分号分隔）和`Cluster`（聚类标签）
   - 脚本自动检测分隔符（逗号、制表符、分号、竖线）

3. **正选择位点列表**（可选但推荐）：在脚本主函数中指定
   - 默认值：`POSITIVE_SELECTION_POSITIONS`（第31-32行的列表）
   - 格式：1-based的氨基酸位置编号列表

### 文件结构示例

```
项目目录/
├── AA_Position_Importance.py            # 主脚本
├── P450_unique_pep_final.fasta          # FASTA序列文件
├── P450_unique_pep_final_esm_clusters_genes.csv  # 聚类标签文件
└── enhanced_position_importance_6repeats/  # 输出目录（运行时创建）
```

## 配置参数

### 脚本内部参数设置

在`main()`函数中（第1089-1097行），可修改以下参数：

```python
# 文件路径配置
FASTA_FILE = "P450_unique_pep_final.fasta"           # FASTA文件路径
CLUSTER_FILE = "P450_unique_pep_final_esm_clusters_genes.csv"  # 聚类文件路径
OUTPUT_DIR = "enhanced_position_importance_6repeats"  # 输出目录名称

# 正选择位点列表（1-based位置编号）
POSITIVE_SELECTION_POSITIONS = [2, 3, 4, 6, 7, 8, 10, 12, 13, 15, 17, 20, 
                                21, 22, 23, 25, 26, 28, 29, 30, 31, 32, 34, 
                                35, 36, 44, 46, 49, 57, 75, 96]

# 分析参数
N_AA = 100  # 分析N端100个氨基酸
```

### ESM2模型配置

在初始化分析器时（第1103-1107行）：

```python
analyzer = EnhancedPositionImportanceAnalyzer(
    model_name="facebook/esm2_t33_650M_UR50D",  # 使用650M参数的ESM2模型
    device="cpu"  # 或 "cuda"（如有GPU）
)
```

### 重复实验设置

- 使用6个不同的随机种子：`[42, 123, 456, 789, 101112, 131415]`
- 每个重复使用5折分层交叉验证
- 关键位置识别标准：平均排名前20%且变异系数<0.5

## 运行脚本

### 基本执行

```bash
# 确保数据文件在同一个目录
python AA_Position_Importance.py
```

### 执行步骤

脚本按以下顺序执行：

1. **步骤1：加载数据** - 读取FASTA和聚类文件，对齐基因列表
2. **步骤2：准备N端序列** - 提取N端前100个氨基酸（或指定长度）
3. **步骤3：运行6次重复分析** - 每个重复进行5折交叉验证
4. **步骤4：与正选择位点比较** - 超几何检验和相似度分析
5. **步骤5：生成可视化图表** - 多种统计图表
6. **步骤6：生成综合报告** - 包含详细分析的文本报告
7. **步骤7：保存详细数据** - CSV和JSON格式的完整结果

### 进度监控

脚本提供详细的进度信息：

- 每个重复的开始和完成状态
- 位置重要性计算进度
- 统计检验结果
- 文件保存确认

## 输出结果

### 目录结构

```
enhanced_position_importance_6repeats/
├── repeat_1_seed_42/                    # 第一次重复结果
│   ├── position_importance.csv          # 单次重复的位置重要性
│   └── metadata.json                    # 单次重复的元数据
├── repeat_2_seed_123/                   # 第二次重复结果
│   ├── position_importance.csv
│   └── metadata.json
├── ...                                  # 重复3-6的结果
├── integrated_results/                  # 整合结果目录
│   ├── integrated_position_importance.csv   # 整合的位置重要性数据
│   ├── comparison_results.json         # 与正选择位点比较结果
│   ├── analysis_metadata.json          # 分析元数据
│   └── plots/                          # 可视化图表目录
│       ├── integrated_position_importance.png/pdf
│       ├── svm_vs_positive_selection.png/pdf
│       └── position_stability.png/pdf
├── comprehensive_analysis_report.txt    # 综合分析报告
└── 控制台输出关键结果和摘要
```

### 关键输出文件说明

1. **整合位置重要性文件** (`integrated_position_importance.csv`)
   - 包含每个位置的详细统计：均值、标准差、中位数、变异系数
   - 标记关键位置和正选择位点
   - 适用于R语言或Python的进一步分析

2. **比较结果文件** (`comparison_results.json`)
   - SVM关键位置与正选择位点的重叠分析
   - 超几何检验结果（p值、期望重叠、实际重叠）
   - Jaccard相似系数和重叠比例

3. **综合分析报告** (`comprehensive_analysis_report.txt`)
   - 8部分结构化报告：执行摘要、关键位置识别、稳定性分析、比较结果、生物学意义、实验建议、方法学评估、结论
   - 包含详细的统计数据和生物学解释

4. **可视化图表** (`plots/`)
   - 整合位置重要性图：条形图带误差条
   - 重复热图：显示6次重复的重要性模式
   - 累积重要性分布：展示最重要的位置贡献
   - SVM与正选择位点比较：重叠可视化
   - 位置稳定性分析：变异系数图

5. **单次重复数据** (`repeat_*/`)
   - 每次独立实验的详细结果
   - 用于检查重复之间的变异性和一致性

## 结果解读指南

### 关键统计指标

1. **位置重要性分数**

   - 范围：0到正无穷，越高表示对分类越重要
   - 解读：分数>平均值+标准差的位置通常具有生物学意义

2. **变异系数 (CV)**

   - 公式：标准差 / 均值
   - 低CV (<0.5)：位置重要性在重复间稳定
   - 高CV (>1.0)：位置重要性不稳定，需谨慎解释

3. **关键位置识别标准**

   ```python
   # 满足以下两个条件的位置被识别为关键位置：
   1. 平均排名在前20%（按重要性降序排列）
   2. 变异系数 < 0.5（在重复间稳定）
   ```

4. **与正选择位点的统计检验**

   - **超几何检验p值**：评估重叠是否显著多于随机

     - p < 0.05：重叠具有统计学显著性
     - p < 0.01：高度显著

   - **Jaccard相似系数**：重叠比例

     - 范围：0（无重叠）到1（完全重叠）

     - >0.2：中度重叠

     - >0.5：高度重叠

### 生物学意义判断

1. **显著重叠** (p < 0.05)
   - 支持假说：功能分类的重要位置也是适应性进化的目标
   - 建议：优先实验验证这些重叠位置

2. **无显著重叠** (p ≥ 0.05)
   - 提示：功能分类和适应性进化可能涉及不同的残基集
   - 可能原因：功能限制、不同的选择压力

3. **位置功能推测**
   - **早期位置 (1-20)**：信号肽、膜靶向、初始折叠
   - **中期位置 (21-40)**：底物通道入口、结构基序
   - **晚期位置 (>40)**：活性位点邻近、辅因子结合

### 实际应用示例

```python
# 示例：识别的高置信度靶点
if p_value < 0.05 and overlap_ratio > 0.3:
    print("✅ 高置信度靶点：功能重要且受正选择")
    # 建议优先进行点突变实验
    
elif key_positions and p_value >= 0.05:
    print("⚠️ 功能重要但与进化压力可能不同")
    # 建议进行功能验证但不排除结构限制
```

## 方法学细节

### ESM2嵌入提取

- 模型：`facebook/esm2_t33_650M_UR50D`（6.5亿参数）
- 嵌入方法：每个位置的独立残基嵌入
- 序列长度：自动截断为1024个token，保留[CLS]和[SEP]标记

### SVM分析与反投影

1. **特征提取**：序列的平均ESM2嵌入
2. **预处理**：标准化 + PCA降维（保留最多50个主成分）
3. **模型训练**：LinearSVM，5折交叉验证，类别平衡权重
4. **权重反投影**：$w_{original} = PCA_{components}^T \cdot w_{PCA}$
5. **位置重要性计算**：$importance_{pos} = mean(|w_{original} \cdot embedding_{pos}|)$

### 稳健性保证

- 6次独立重复，不同随机种子
- 5折交叉验证，分层采样
- PCA降维提高数值稳定性
- 变异系数过滤确保结果稳定

## 故障排除

### 常见问题与解决方案

1. **内存不足错误**

```python
# 解决方案1：使用较小的ESM2模型
analyzer = EnhancedPositionImportanceAnalyzer(
    model_name="facebook/esm2_t12_35M_UR50D",  # 3500万参数
    device="cpu"
)

# 解决方案2：减小批处理大小（修改第125行）
batch_size = 1  # 原为2
```

2. **JSON序列化错误**（已修复）

- 脚本已内置`NumpyEncoder`和`convert_to_serializable`函数
- 自动处理NumPy数据类型、NaN值和复杂对象
- 如果仍出现问题，检查Python的json模块版本

3. **数据加载失败**

```bash
# 检查文件格式
# FASTA：确保标准格式，无空行或特殊字符
# CSV：确保Gene_List和Cluster列存在
# 基因名匹配：确保FASTA头与CSV中的基因标识符一致
```

4. **正选择位点列表修改**

```python
# 在main()函数中修改POSITIVE_SELECTION_POSITIONS
# 使用1-based位置编号，如[1, 2, 3, 5, 8, 13]
POSITIVE_SELECTION_POSITIONS = [您的位点列表]
```

5. **可视化错误**

```bash
# 安装缺失的包
pip install seaborn
pip install matplotlib

# 或禁用特定可视化（注释相关代码段）
```

### 错误信息与解决方案

1. **"No common genes found"**
   - 原因：FASTA和CSV文件的基因标识符不匹配
   - 解决：统一基因命名格式，或修改解析逻辑

2. **"CUDA out of memory"**
   - 原因：GPU内存不足
   - 解决：使用CPU版本或减小批处理大小

3. **"ValueError: shapes not aligned"**
   - 原因：特征维度不匹配
   - 解决：检查PCA组件数与特征维度，或重启脚本

## 高级定制

### 修改分析参数

```python
# 1. 修改关键位置识别标准（第643-646行）
top_percentage = 0.15  # 原为0.2（前20%）
cv_threshold = 0.4     # 原为0.5

# 2. 修改重复次数和种子（第407-408行）
random_seeds = [42, 84, 168, 336, 672, 1344]  # 自定义种子

# 3. 修改SVM参数（第204-210行）
svm = LinearSVC(
    C=0.1,              # 正则化强度（原为1.0）
    class_weight='balanced',
    max_iter=5000,      # 最大迭代次数（原为10000）
    random_state=random_state,
    dual=True,          # 对偶问题（原为False）
    tol=1e-3            # 容忍度（原为1e-4）
)
```

### 扩展分析功能

```python
# 添加新的可视化类型
def custom_visualization(self, integrated_results, output_dir):
    # 自定义绘图代码
    pass

# 在main()函数中调用
analyzer.custom_visualization = custom_visualization.__get__(analyzer)
```

### 批量分析多个区域

```python
# 修改main()函数以分析多个长度
for n_aa in [50, 75, 100, 125]:
    print(f"分析N端 {n_aa}AA")
    analyzer.six_repeats_analysis(sequences, labels, n_aa=n_aa)
```

## 科学验证与假设检验

### 统计验证方法

1. **重复间一致性**：计算6次重复的相关系数矩阵
2. **排名稳定性**：Kendall's W一致性系数
3. **与进化数据整合**：超几何检验评估功能-进化关联
4. **效应量评估**：Jaccard相似系数和重叠比例

### 假说检验框架

```python
# 零假设H0：SVM关键位置与正选择位点的重叠是随机的
# 备择假设H1：重叠显著多于随机（p < 0.05）

if p_value < 0.05:
    print("拒绝H0：功能重要位置与适应性进化位置显著重叠")
    print("支持P450功能-进化协同进化的假说")
else:
    print("无法拒绝H0：重叠不显著")
    print("功能分类和适应性进化可能涉及不同的机制")
```

## 引用与致谢

如使用本脚本，请引用：

- **ESM2模型**：Lin et al. (2022) *Nature Methods*
- **SVM方法**：Cortes & Vapnik (1995) *Machine Learning*
- **正选择分析**：Yang (2007) *Molecular Biology and Evolution*

**方法学参考文献**：  
*基于重复线性SVM和ESM2嵌入的植物P450蛋白位置重要性分析：功能-进化整合框架*

## 技术支持与联系

技术问题解决步骤：

1. 检查控制台错误信息和堆栈跟踪
2. 验证数据文件格式和完整性
3. 检查依赖包版本兼容性
4. 尝试简化测试数据集

科学研究问题：

- 结果解读：参考`comprehensive_analysis_report.txt`中的生物学解释
- 方法选择：根据研究目标和数据特性调整参数
- 实验设计：基于关键位置列表设计验证实验

脚本更新日志：

- v1.0：初始版本，6次重复分析框架
- v1.1：修复JSON序列化问题，增强稳定性
- v1.2：添加综合报告生成，改进可视化

---

**脚本版本**：AA_Position_Importance.py v1.1  
**最后更新**：2025年  
**主要特点**：6次重复稳健分析、正选择位点整合、修复JSON序列化  
**应用场景**：植物P450蛋白功能位点识别、功能-进化关联研究、实验靶点优先级排序

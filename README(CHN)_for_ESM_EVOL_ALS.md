# ESM Evolutionary Analysis Tool (ESM进化分析工具)
**作者:** **Qiao Weizhuo**  
**许可License:** **MIT**

## 概述 | Overview

这是一个基于ESM2（Evolutionary Scale Modeling 2）蛋白质语言模型的进化分析工具，用于比较蛋白质序列的功能聚类与系统发育聚类之间的差异。该工具能够识别在功能上趋同进化但系统发育关系较远的蛋白质序列。

## 主要功能 | Key Features

### 1. 特征提取 | Feature Extraction
- 使用ESM2模型（650M或15B参数）提取蛋白质序列的深度表征
- 支持批量处理，高效提取嵌入向量
- 自动保存特征向量为CSV格式

### 2. 聚类分析 | Clustering Analysis
- **层次聚类**：支持多种连接方法（average, ward, complete, single）和距离度量（cosine, euclidean, correlation）
- **系统发育树聚类**：基于Newick格式的系统发育树进行单系群聚类
- **动态规划聚类**：优化系统发育树的聚类分割
- **非单系群修复**：自动检测和修复非单系群簇

### 3. 比较分析 | Comparative Analysis
- 计算V-measure和调整互信息(AMI)评估聚类一致性
- 识别功能与系统发育不一致的"分歧序列"
- 比较ESM距离矩阵与系统发育距离矩阵

### 4. 降维可视化 | Dimensionality Reduction and Visualization
- UMAP降维（2D）
- PCoA（主坐标分析）比较
- 多种可视化图表生成

### 5. 数据导出 | Data Export
- 桑基图数据表（用于ESM与系统发育聚类关系可视化）
- 距离矩阵（CSV格式）
- 聚类分配结果
- Newick格式的层次聚类树

## 文件结构 | File Structure

```
esm_evolution_tool/
├── esm_evol_als.py           # 主程序入口
├── core/                     # 核心功能模块
│   ├── __init__.py
│   ├── data_processor.py    # 数据预处理和特征提取
│   ├── clustering.py        # 聚类算法实现
│   ├── metrics.py           # 评估指标计算
│   └── utils.py             # 工具函数
└── visualization/           # 可视化模块
    ├── __init__.py
    ├── exporter.py          # 数据导出功能
    └── plotter.py           # 图表绘制功能
```

## 安装依赖 | Installation Requirements

```bash
# 基础依赖
pip install numpy scipy pandas matplotlib seaborn scikit-learn

# ESM2模型相关
pip install torch transformers

# 生物信息学工具
pip install biopython ete3

# 降维和可视化
pip install umap-learn plotly

# 可选：用于树可视化的额外包
pip install cairosvg  # 用于系统发育树导出为PNG
```

## 使用方法 | Usage

### 基本命令 | Basic Command

```bash
python esm_evol_als.py -i input.fasta -t tree.nwk
```

### 完整参数 | Complete Parameters

```bash
python esm_evol_als.py \
  -i input.fasta \               # 输入FASTA文件（必需）
  -t tree.nwk \                  # 系统发育树文件（必需，Newick格式）
  -m esm2-650M \                 # ESM2模型选择：esm2-650M 或 esm2-15B
  -l average \                   # 层次聚类连接方法
  --metric cosine \              # 距离度量方法
  -d 0.5 \                       # 层次聚类距离阈值（可选）
  -p 10                          # 系统发育树聚类数量（可选）
```

### 参数说明 | Parameter Description

- `-i, --input`: 输入的蛋白质序列FASTA文件
- `-t, --tree`: Newick格式的系统发育树文件
- `-m, --model`: ESM2模型选择（默认：esm2-650M）
  - `esm2-650M`: 650M参数模型，速度较快
  - `esm2-15B`: 15B参数模型，精度更高但需要更多资源
- `-l, --linkage`: 层次聚类连接方法
  - `average`: 平均连接（默认）
  - `ward`: Ward方法
  - `complete`: 完全连接
  - `single`: 单连接
- `--metric`: 距离度量方法
  - `cosine`: 余弦距离（默认）
  - `euclidean`: 欧氏距离
  - `correlation`: 相关系数距离
- `-d, --distance-threshold`: 层次聚类距离阈值，用于确定聚类数量
- `-p, --phylo-clusters`: 系统发育树聚类数量，如不指定则使用与ESM聚类相同的数量

## 输出文件 | Output Files

运行程序后会生成以下文件（以输入文件名为前缀）：

### 数据文件 | Data Files
1. `{prefix}_esm2_embeddings.csv` - ESM2嵌入特征矩阵
2. `{prefix}_cluster_assignments.csv` - 聚类分配结果
3. `{prefix}_distance_matrix.csv` - ESM距离矩阵
4. `{prefix}_phylogenetic_distance_matrix.csv` - 系统发育距离矩阵
5. `{prefix}_evolution_analysis.csv` - 进化分析结果
6. `{prefix}_divergent_sequences.csv` - 分歧序列列表
7. `{prefix}_pcoa_results.csv` - PCoA分析结果
8. `{prefix}_linkage_matrix.csv` - 层次聚类连接矩阵

### 桑基图数据 | Sankey Diagram Data
1. `{prefix}_esm_clusters_genes.csv` - ESM聚类基因列表
2. `{prefix}_phylogenetic_clusters_genes.csv` - 系统发育聚类基因列表
3. `{prefix}_sankey_flow_genes.csv` - 桑基图流量基因列表
4. `{prefix}_sankey_detailed_assignments.csv` - 详细分配数据

### 可视化图表 | Visualization Charts
1. `{prefix}_hclust_dendrogram.png` - 层次聚类树状图
2. `{prefix}_phylogeny_tree.png` - 系统发育树图
3. `{prefix}_cluster_comparison.png` - 聚类比较图
4. `{prefix}_convergence_heatmap.png` - 趋同进化热图
5. `{prefix}_pcoa_comparison.png` - PCoA比较图
6. `{prefix}_distance_heatmap.png` - 距离矩阵热图
7. `{prefix}_hclust_tree.nwk` - Newick格式的层次聚类树
8. `{prefix}_sankey_diagram.html` - 交互式桑基图（HTML格式）

## 算法细节 | Algorithm Details

### 1. 特征提取 | Feature Extraction
- 使用ESM2模型的最后一层CLS token作为序列表征
- 自动处理序列长度差异（padding和truncation）
- 支持GPU加速（如可用）

### 2. 层次聚类 | Hierarchical Clustering
- 基于scipy的层次聚类实现
- 支持多种距离度量和连接方法
- 自动确定聚类阈值（70%最大距离）

### 3. 系统发育树聚类 | Phylogenetic Tree Clustering
- 动态规划算法寻找最优单系群分割
- 非单系群检测和自动修复
- 确保所有簇都是单系群

### 4. 差异分析 | Divergence Analysis
- 使用V-measure和调整互信息(AMI)评估聚类一致性
- 识别功能与系统发育不一致的序列
- 提供详细的差异序列列表

## 应用场景 | Application Scenarios

1. **趋同进化研究**：识别功能相似但系统发育关系较远的蛋白质
2. **功能注释验证**：验证基于序列相似性的功能注释准确性
3. **蛋白质工程**：识别具有相似功能但不同进化起源的蛋白质变体
4. **比较基因组学**：研究不同物种间蛋白质功能的保守性和分化

## 注意事项 | Notes

1. **硬件要求**：
   - esm2-650M模型需要约2.5GB GPU内存
   - esm2-15B模型需要约30GB GPU内存
   - 大型数据集可能需要大量内存和计算时间

2. **输入数据要求**：
   - FASTA文件中的序列ID应与系统发育树中的叶节点名称一致
   - 系统发育树应为有根树

3. **可视化依赖**：
   - 桑基图需要plotly库，如未安装将跳过生成
   - 系统发育树可视化需要ete3的树渲染功能

4. **性能优化**：
   - 对于大型数据集，建议使用较小的模型或增加批处理大小
   - 层次聚类的时间复杂度为O(n²)，大型数据集可能较慢

## 引用 | Citation

如使用本工具，请引用相关依赖库：
- ESM2: [Evolutionary Scale Modeling](https://github.com/facebookresearch/esm)
- scikit-learn: [Machine Learning in Python](https://scikit-learn.org/)
- ete3: [Python Environment for Tree Exploration](http://etetoolkit.org/)

## 许可证 | License

此工具基于开源软件构建，请遵循各依赖库的相应许可证。

## 支持与贡献 | Support and Contribution

如遇到问题或有改进建议，请提交issue或pull request。

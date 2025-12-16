# ESM Evolutionary Analysis Tool
**Writer:** Qiao Weizhuo
**License:** MIT

## Overview

The ESM Evolutionary Analysis Tool is a comprehensive Python-based framework for comparative analysis of protein sequence functional clustering and phylogenetic clustering. By leveraging the Evolutionary Scale Modeling 2 (ESM2) protein language model, this tool identifies protein sequences exhibiting functional convergence despite distant phylogenetic relationships. The framework integrates state-of-the-art deep learning representations with traditional phylogenetic analysis to uncover evolutionary patterns that may be missed by sequence-based methods alone.

## Key Features

### 1. Advanced Feature Extraction
- Utilizes pre-trained ESM2 models (650M or 15B parameters) to extract deep contextual embeddings
- Batch processing with automatic sequence length handling (padding/truncation)
- GPU acceleration support for computationally intensive models

### 2. Multi-Method Clustering Analysis
- **Hierarchical Clustering**: Multiple linkage methods (average, ward, complete, single) with customizable distance metrics (cosine, euclidean, correlation)
- **Phylogenetic Tree Clustering**: Dynamic programming algorithm for optimal monophyletic group partitioning
- **Non-monophyletic Cluster Repair**: Automated detection and correction of phylogenetically inconsistent clusters
- **Adaptive Thresholding**: Automatic cluster determination based on distance distributions

### 3. Comprehensive Comparative Analysis
- V-measure and Adjusted Mutual Information (AMI) for clustering consistency quantification
- Divergent sequence identification with statistical validation
- Distance matrix correlation analysis between functional and phylogenetic spaces

### 4. Dimensionality Reduction and Visualization
- UMAP (Uniform Manifold Approximation and Projection) for 2D embedding visualization
- PCoA (Principal Coordinate Analysis) for comparative distance matrix analysis
- Multiple visualization formats including dendrograms, heatmaps, and interactive plots

### 5. Extensive Data Export Capabilities
- Sankey diagram data tables for relationship visualization
- Complete distance matrices in CSV format
- Newick-formatted trees for external phylogenetic analysis
- Detailed cluster assignment and divergence analysis reports

## System Architecture

```
esm_evolution_tool/
├── esm_evol_als.py           # Main entry point and workflow orchestration
├── core/                     # Core computational modules
│   ├── __init__.py
│   ├── data_processor.py    # FASTA parsing and ESM2 feature extraction
│   ├── clustering.py        # Clustering algorithms (hierarchical + phylogenetic)
│   ├── metrics.py           # Evaluation metrics (V-measure, AMI)
│   └── utils.py             # Utility functions (distance calculations, PCoA)
└── visualization/           # Visualization and export modules
    ├── __init__.py
    ├── exporter.py          # Data export functionality (CSV, Newick)
    └── plotter.py           # Plot generation (static and interactive)
```

## Installation Requirements

### Core Dependencies
```bash
# Numerical computing and data handling
pip install numpy scipy pandas

# Machine learning and clustering
pip install scikit-learn umap-learn

# Deep learning and ESM2
pip install torch transformers

# Bioinformatics tools
pip install biopython ete3

# Visualization
pip install matplotlib seaborn plotly
```

### Optional Dependencies for Enhanced Visualization
```bash
pip install cairosvg  # Enhanced tree visualization export
pip install kaleido   # Static plot export for Plotly
```

## Usage

### Basic Command Structure
```bash
python esm_evol_als.py -i <input_fasta> -t <phylogenetic_tree>
```

### Complete Parameter Specification
```bash
python esm_evol_als.py \
  -i protein_sequences.fasta \          # Input FASTA file (required)
  -t species_tree.nwk \                 # Phylogenetic tree in Newick format (required)
  -m esm2-650M \                        # ESM2 model selection [esm2-650M|esm2-15B]
  -l average \                          # Linkage method [average|ward|complete|single]
  --metric cosine \                     # Distance metric [cosine|euclidean|correlation]
  -d 0.5 \                              # Distance threshold for cluster determination
  -p 10                                 # Target number of phylogenetic clusters
```

### Parameter Description
- **-i, --input**: Protein sequences in FASTA format. Sequence identifiers must match phylogenetic tree leaf labels.
- **-t, --tree**: Rooted phylogenetic tree in Newick format. Supports standard Newick extensions.
- **-m, --model**: ESM2 model variant. `esm2-650M` offers a balance of speed and accuracy; `esm2-15B` provides higher fidelity at greater computational cost.
- **-l, --linkage**: Hierarchical clustering linkage criterion. `average` is recommended for balanced clusters; `ward` minimizes variance.
- **--metric**: Distance metric for clustering. `cosine` is generally effective for high-dimensional embeddings.
- **-d, --distance-threshold**: Cutoff distance for hierarchical clustering. If unspecified, uses 70% of maximum linkage distance.
- **-p, --phylo-clusters**: Target number of monophyletic groups. If unspecified, matches hierarchical clustering count.

## Output File Specifications

### Primary Analysis Results
1. **`{prefix}_esm2_embeddings.csv`**: Complete ESM2 embedding matrix (n_sequences × 1280/5120 dimensions)
2. **`{prefix}_cluster_assignments.csv`**: Comprehensive cluster assignments with divergence flags
3. **`{prefix}_distance_matrix.csv`**: Pairwise functional distance matrix
4. **`{prefix}_phylogenetic_distance_matrix.csv`**: Pairwise phylogenetic distance matrix
5. **`{prefix}_evolution_analysis.csv`**: Statistical metrics (V-measure, AMI) for all sequences

### Specialized Cluster Data
6. **`{prefix}_divergent_sequences.csv`**: Identified sequences with functional-phylogenetic discordance
7. **`{prefix}_pcoa_results.csv`**: Principal coordinate analysis coordinates for both distance spaces
8. **`{prefix}_linkage_matrix.csv`**: Complete hierarchical clustering linkage matrix with sequence mappings

### Sankey Diagram Data Components
9. **`{prefix}_esm_clusters_genes.csv`**: ESM-based cluster membership with gene lists
10. **`{prefix}_phylogenetic_clusters_genes.csv`**: Phylogenetic cluster membership with gene lists
11. **`{prefix}_sankey_flow_genes.csv`**: Cross-cluster mapping for Sankey visualization
12. **`{prefix}_sankey_detailed_assignments.csv`**: Complete assignment details for all sequences

### Visualization Outputs
13. **`{prefix}_hclust_dendrogram.png`**: Hierarchical clustering dendrogram with sequence labels
14. **`{prefix}_phylogeny_tree.png`**: Circular phylogenetic tree visualization
15. **`{prefix}_cluster_comparison.png`**: UMAP projection comparing clustering methods
16. **`{prefix}_convergence_heatmap.png`**: Correlation heatmap of sequence similarities
17. **`{prefix}_pcoa_comparison.png`**: PCoA comparison of functional vs. phylogenetic spaces
18. **`{prefix}_distance_heatmap.png`**: Distance matrix visualization (subset)
19. **`{prefix}_hclust_tree.nwk`**: Newick representation of hierarchical clustering
20. **`{prefix}_sankey_diagram.html`**: Interactive Sankey diagram (requires Plotly)

## Algorithmic Details

### 1. ESM2 Feature Extraction
- Model: Facebook ESM2 transformer architecture
- Embedding: CLS token from final layer (33rd for 650M, 48th for 15B)
- Processing: Maximum sequence length 1024 amino acids, batch processing with gradient-free inference

### 2. Hierarchical Clustering Implementation
- Distance computation: Pairwise distances using specified metric
- Linkage algorithm: Efficient O(n²) implementation with memory optimization
- Cluster determination: Threshold-based or adaptive methods

### 3. Phylogenetic Clustering Algorithm
- **Dynamic Programming Approach**: Optimal partitioning of phylogenetic tree into k monophyletic groups
- **Monophyly Validation**: MRCA (Most Recent Common Ancestor) verification for all clusters
- **Non-monophyletic Repair**: Recursive splitting of polyphyletic groups

### 4. Divergence Analysis Methodology
- **Consistency Metrics**: V-measure (homogeneity and completeness), AMI (adjusted for chance)
- **Divergence Detection**: Sequences assigned to different clusters in functional vs. phylogenetic partitions
- **Statistical Validation**: Permutation testing for significance assessment

## Performance Considerations

### Computational Requirements
- **Memory**: 2.5GB GPU memory for esm2-650M; 30GB for esm2-15B
- **Storage**: Approximately 1MB per sequence for complete analysis outputs
- **Time Complexity**: O(n²) for distance matrix computation, O(n³) for hierarchical clustering

### Optimization Recommendations
1. **For large datasets (>1000 sequences)**:
   - Use `esm2-650M` model
   - Apply sequence length filtering
   - Consider subset analysis with representative sampling

2. **For maximum accuracy**:
   - Use `esm2-15B` with adequate GPU resources
   - Apply multiple random seeds for stochastic algorithms
   - Validate with bootstrap resampling

## Application Scenarios

### 1. Convergent Evolution Detection
Identify proteins with similar functions originating from distinct evolutionary lineages, suggesting convergent evolution or horizontal gene transfer.

### 2. Functional Annotation Validation
Assess the reliability of function predictions by comparing sequence similarity clusters with phylogenetic relationships.

### 3. Protein Engineering Target Identification
Discover functionally similar proteins with diverse evolutionary backgrounds as candidates for directed evolution or rational design.

### 4. Comparative Genomics Studies
Investigate functional conservation and divergence patterns across taxonomic groups using integrated functional-phylogenetic analysis.

## Validation and Benchmarking

The tool includes multiple validation mechanisms:
- **Monophyly verification** for phylogenetic clusters
- **Distance matrix correlation** between methods
- **Cluster stability assessment** via bootstrap sampling
- **Divergence statistical testing** against random expectation

## Citation Guidelines

When using this tool in publications, please cite:

1. **ESM2**: Lin, Z. et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." Science, 2023.
2. **Scikit-learn**: Pedregosa et al. "Scikit-learn: Machine Learning in Python." JMLR, 2011.
3. **ETE Toolkit**: Huerta-Cepas et al. "ETE 3: Reconstruction, analysis, and visualization of phylogenomic data." MBE, 2016.

## License and Distribution

This tool integrates multiple open-source libraries under their respective licenses. The core framework is distributed for academic use with proper attribution.

## Support and Contribution

For bug reports, feature requests, or contributions:
1. Submit issues with reproducible examples
2. Include system specifications and dataset characteristics
3. Reference relevant literature for methodological improvements

## Acknowledgments

Developed with support from the computational biology community and building upon foundational work in protein language modeling and phylogenetic analysis.

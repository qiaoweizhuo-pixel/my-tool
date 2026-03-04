**Classifier_Selection.py - ReadMe (English Version)**

## Script Function Overview

This script is designed to systematically compare the performance differences of various sequence regions (N-terminal, C-terminal, middle region, random fragments, full-length sequence) of plant P450 proteins in functional classification tasks. By integrating the ESM2 protein language model, multiple machine learning classifiers, and dimensionality reduction visualization methods, this script aims to validate the specificity of N-terminal sequence signals, distinguish the contributions of sequence order information from amino acid composition information, and provide quantitative evidence for understanding the structural basis of P450 protein function evolution.

## Main Functional Modules

1. **Multi-Region Feature Extraction**: Extracts ESM2 embedding vectors for N-terminal, C-terminal, middle region, random fragments, and full-length sequences.
2. **Control Experiment Design**:
   - Shuffled N-terminal: Shuffles N-terminal amino acid order to control for amino acid composition.
   - Random Fragments: Extracts fragments from random positions within the sequence.
   - Multiple Random Fragment Trials: Assesses the robustness of random position selection.
3. **Classifier Comparison**: Compares the classification performance of LinearSVM, LogisticRegression, and RandomForest.
4. **Multi-Dimensional Visualization**: Provides results from three dimensionality reduction methods: PCA, t-SNE, and UMAP.
5. **Statistical Validation**: Includes bootstrap confidence intervals, paired t-tests, and effect size analysis.
6. **Data Export**: Exports all raw data required for visualization, facilitating secondary analysis in R.

## System Requirements and Dependencies

### Python Version Requirements
- Python 3.8 or higher

### Required Python Packages
```bash
pip install torch transformers numpy pandas matplotlib seaborn scikit-learn tqdm scipy
```

### Optional Python Packages (for advanced dimensionality reduction)
```bash
# For t-SNE (already included in scikit-learn)
# For UMAP (install if needed)
pip install umap-learn
```

## File Preparation

### Required Data Files
Place the following files in the same directory as the script:

1. **FASTA File**: Contains P450 protein sequences.
   - Default filename: `P450_unique_pep_final.fasta`
   - Format: Standard FASTA format, with headers containing gene identifiers.

2. **Cluster Label File**: Contains the mapping from genes to functional categories/subfamilies.
   - Default filename: `P450_unique_pep_final_esm_clusters_genes.csv`
   - Required columns: `Gene_List` (list of gene identifiers, separated by semicolons) and `Cluster` (category labels).

### Example File Structure
```
Project Directory/
├── Classifier_Selection.py          # Main script
├── P450_unique_pep_final.fasta      # FASTA sequence file
├── P450_unique_pep_final_esm_clusters_genes.csv  # Cluster label file
└── enhanced_p450_region_analysis/   # Output directory generated after running the script
```

## Configuration Parameters

### Key Parameter Descriptions (located in the `main()` function)
```python
# File path configuration
FASTA_FILE = "P450_unique_pep_final.fasta"
CLUSTER_FILE = "P450_unique_pep_final_esm_clusters_genes.csv"
OUTPUT_DIR = "enhanced_p450_region_analysis"

# Analyzer initialization parameters
analyzer = EnhancedP450Analyzer(
    model_name="facebook/esm2_t33_650M_UR50D",  # ESM2 model version
    n_aa=100,      # Length of N-terminal amino acids to analyze
    c_aa=100,      # Length of C-terminal amino acids to analyze
    middle_aa=100, # Length of the middle region to analyze
    random_aa=100, # Length of random fragments
    device="cpu",  # Use CPU; change to "cuda" if GPU is available
    use_tsne=True, # Enable t-SNE dimensionality reduction
    use_umap=True, # Enable UMAP dimensionality reduction
    num_random_trials=20  # Number of random fragment trials
)
```

## Running the Script

1. **Ensure Data Files Are in Place**: Place the FASTA file and cluster file in the same directory as the script.

2. **Run the Analysis**:
```bash
python Classifier_Selection.py
```

3. **Monitor Output**: Progress information and intermediate results will be displayed during script execution.

## Output Results

### Directory Structure
```
enhanced_p450_region_analysis/
├── plots/                                      # Visualization plot directory
│   ├── classifier_performance_comparison_bar.png/pdf
│   ├── region_performance_comparison_bar.png/pdf
│   ├── performance_heatmap.png/pdf
│   ├── n_terminal_vs_shuffled_comparison.png/pdf
│   ├── region_dimension_reduction.png/pdf
│   └── random_fragment_distributions.png/pdf
├── *.csv                                      # Raw data files
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
└── Console output of statistical results and key conclusions
```

### Main Output File Descriptions

1. **Projection Data Files** (`*_projection_data.csv`)
   - Contain 2D coordinates (PCA/t-SNE/UMAP) for each region.
   - Used for generating advanced custom visualizations in R.

2. **Classifier Performance File** (`classifier_performance.csv`)
   - Cross-validation accuracy statistics for each region and each classifier.

3. **Detailed Accuracy Data** (`detailed_accuracy_data.csv`)
   - Accuracy for each cross-validation fold, suitable for variance analysis.

4. **Region Performance Comparison** (`region_performance_comparison.csv`)
   - Direct comparison of each region under the same cross-validation splits.

5. **Statistical Summary File** (`summary_statistics.csv`)
   - Summarized statistical data for generating bar charts.

## Results Interpretation Guide

### Key Statistical Metrics

1. **Comparison of N-terminal vs. Shuffled N-terminal**
   - N-terminal > Shuffled N-terminal: Indicates the presence of sequence order information.
   - N-terminal ≈ Shuffled N-terminal: Only amino acid composition information is at play.

2. **Comparison of N-terminal vs. Other Regions**
   - N-terminal performs best: Supports the hypothesis of N-terminal signal specificity.
   - Multiple regions perform similarly: Functional information is widely distributed.

3. **Multiple Random Fragment Trials**
   - Stable performance of random fragments: Indicates robustness of the classification task.
   - High variability in random fragment performance: Results may be influenced by random position selection.

### Biological Significance Judgment Criteria

1. **Strong Specificity Signal**:
   ```
   N-terminal > Shuffled N-terminal
   N-terminal > C-terminal
   N-terminal > Middle Region
   N-terminal > Random Fragments
   ```

2. **Sequence Order Dependence**:
   ```
   N-terminal > Shuffled N-terminal
   Cohen's d > 0.5 (moderate or higher effect size)
   p-value < 0.05 (statistically significant)
   ```

## Troubleshooting

### Common Issues and Solutions

1. **Insufficient Memory**
   - Reduce batch size: Modify the `batch_size` parameter.
   - Use a smaller ESM2 model: e.g., `facebook/esm2_t12_35M_UR50D`.
   - Enable PCA dimensionality reduction: Reduce feature dimensions.

2. **Data Loading Failure**
   - Check FASTA file format: Ensure it is standard FASTA format.
   - Check CSV delimiter: The script detects automatically, but can be manually specified.
   - Verify gene identifier consistency: Ensure gene names match between the FASTA and CSV files.

3. **Slow or Failed Model Download**
   - Download the model manually: Specify a local path after downloading from Hugging Face.
   - Use a proxy: Set up a network proxy.
   - Choose a smaller model: e.g., `facebook/esm2_t6_8M_UR50D`.

4. **UMAP/t-SNE Unavailable**
   - Install missing packages: `pip install umap-learn`.
   - Disable corresponding features: Set `use_umap=False` or `use_tsne=False`.

## Advanced Customization

### Modifying Analyzed Region Lengths
```python
# Modify in the main() function
analyzer = EnhancedP450Analyzer(
    n_aa=50,      # Analyze the first 50 N-terminal amino acids
    c_aa=50,      # Analyze the last 50 C-terminal amino acids
    # ... other parameters
)
```

### Adding Classifiers
```python
# Add classifiers in the compare_feature_regions method
classifiers = {
    'LinearSVM': LinearSVC(...),
    'LogisticRegression': LogisticRegression(...),
    'RandomForest': RandomForest(...),
    'New_Classifier': Classifier_Instance
}
```

### Adjusting Statistical Thresholds
```python
# Modify statistical significance thresholds
alpha = 0.01  # Stricter significance level
# Modify in the paired t-test section
```

## Citations and Acknowledgments

If using this script, please cite the relevant methods:
- ESM2 protein language model: Lin et al. (2022) Nature Methods
- Machine learning framework: scikit-learn
- Dimensionality reduction methods: PCA (scikit-learn), t-SNE (van der Maaten & Hinton, 2008), UMAP (McInnes et al., 2018)

## Technical Support and Feedback

For issues or suggestions, please:
1. Check console error messages.
2. Review log files in the output directory.
3. Ensure all dependency package versions are compatible.
4. Provide reproducible example data for debugging.

---
*Script Designer: P450 Evolution Analysis Project Team*  
*Last Updated: 2025*  
*Version: EnhancedP450Analyzer v1.0*
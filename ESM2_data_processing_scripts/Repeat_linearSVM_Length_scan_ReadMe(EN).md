# Repeat_linearSVM_Length_scan.py - Plant P450 N-terminal Length Scanning Analysis with Repeated Linear SVM

## Script Function Overview

This script performs a comprehensive analysis of N-terminal length effects on plant P450 protein classification performance. Using a repeated Linear SVM approach (6 repeats by default), it systematically scans N-terminal lengths from 40 to 100 amino acids to validate three evolutionary hypotheses: rapid divergence (40-49AA), transition (50-84AA), and noise removal (85-100AA). The script integrates ESM2 embeddings with robust machine learning to ensure statistically reliable results.

## Main Functional Modules

1. **Repeated Analysis Framework**: Executes 6 independent Linear SVM analyses with different random seeds to ensure result robustness
2. **Multi-length Feature Extraction**: Extracts ESM2 embeddings for N-terminal sequences of varying lengths (40-100AA)
3. **Hypothesis Testing**: Specifically tests three evolutionary periods in P450 N-terminal evolution
4. **Statistical Integration**: Aggregates results across 6 independent repeats to calculate mean performance, significance ratios, and effect sizes
5. **Comprehensive Visualization**: Generates integrated plots showing performance curves, significance heatmaps, and turning point analysis
6. **Data Export**: Exports detailed results in CSV and JSON formats for further analysis

## System Requirements and Dependencies

### Python Version Requirements
- Python 3.8 or higher

### Required Python Packages
```bash
pip install torch transformers numpy pandas matplotlib seaborn scikit-learn tqdm scipy
```

### For Advanced Users
```bash
# Optional: For GPU acceleration (requires CUDA-compatible GPU)
# Ensure appropriate CUDA version is installed
pip install torch --index-url https://download.pytorch.org/whl/cu118  # Example for CUDA 11.8
```

## File Preparation

### Required Data Files
Place the following files in the same directory as the script:

1. **FASTA File**: Contains P450 protein sequences
   - Default filename: `P450_unique_pep_final.fasta`
   - Format: Standard FASTA format with gene identifiers in headers

2. **Cluster Label File**: Contains functional clustering information
   - Default filename: `P450_unique_pep_final_esm_clusters_genes.csv`
   - Required columns: Must include gene identifiers and cluster labels
   - The script automatically detects separators (comma, tab, semicolon, pipe)

### File Structure Example
```
Project Directory/
├── repeat_linearSVM_Length_scan.py      # Main script
├── P450_unique_pep_final.fasta          # FASTA sequence file
├── P450_unique_pep_final_esm_clusters_genes.csv  # Cluster label file
└── p450_length_scan_repeated_svm/       # Output directory (created during runtime)
```

## Configuration Parameters

### Command Line Arguments
The script supports the following command line arguments:

```python
# Basic usage with default parameters
python repeat_linearSVM_Length_scan.py

# Customized analysis
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

### Parameter Descriptions
```python
# Default values and descriptions:
--fasta: "P450_unique_pep_final.fasta"           # FASTA file path
--cluster: "P450_unique_pep_final_esm_clusters_genes.csv"  # Cluster file path
--output: "p450_length_scan_repeated_svm"        # Output directory name
--min_length: 40                                 # Minimum N-terminal length to analyze
--max_length: 100                                # Maximum N-terminal length to analyze
--n_repeats: 6                                   # Number of independent repeats
--device: "cpu"                                  # Compute device ("cpu" or "cuda")
--seed_start: 42                                 # Starting random seed (increments by 100)
```

## Running the Script

### Basic Execution
```bash
# 1. Ensure data files are in the same directory
# 2. Run the script with default parameters
python repeat_linearSVM_Length_scan.py
```

### Advanced Execution
```bash
# Run with custom parameters
python repeat_linearSVM_Length_scan.py \
    --min_length 30 \
    --max_length 150 \
    --n_repeats 10 \
    --device "cuda"  # Use GPU if available
```

### Monitoring Progress
The script provides detailed progress information:
- Model loading status
- Data loading statistics
- Individual repeat progress (1-6)
- Statistical integration steps
- Visualization generation

## Output Results

### Directory Structure
```
p450_length_scan_repeated_svm/                # Main output directory
├── repeat_1_seed_42/                         # Individual repeat results
│   ├── n_terminal_length_scan.csv
│   └── full_length_performance.csv
├── repeat_2_seed_142/                        # Second repeat
│   ├── n_terminal_length_scan.csv
│   └── full_length_performance.csv
├── ...                                       # More repeats (3-6)
├── integrated_results/                       # Aggregated results
│   ├── integrated_n_terminal_length_scan.csv
│   ├── integrated_full_length_performance.csv
│   ├── detailed_all_repeats_accuracy.csv
│   ├── all_results_summary.json
│   └── plots/                                # Visualization directory
│       ├── integrated_length_scan_analysis.png/pdf
│       ├── hypothesis_testing_turning_point.png/pdf
│       └── individual_repeats_comparison.png/pdf
└── console_output.txt                        # Runtime log (if redirected)
```

### Key Output Files

1. **Integrated Analysis Files** (`integrated_results/`)
   - `integrated_n_terminal_length_scan.csv`: Mean and variation across repeats for each length
   - `detailed_all_repeats_accuracy.csv`: Raw accuracy values for each repeat and length
   - `all_results_summary.json`: Comprehensive summary including turning point analysis

2. **Visualization Files** (`plots/`)
   - `integrated_length_scan_analysis.png/pdf`: Four-panel figure showing:
     - Average accuracy curve with standard deviation
     - Significance ratio heatmap
     - Performance gap percentage
     - Accuracy distribution boxplots
   - `hypothesis_testing_turning_point.png/pdf`: Main hypothesis testing figure
   - `individual_repeats_comparison.png/pdf`: Comparison of all 6 repeats

3. **Individual Repeat Files** (`repeat_*/`)
   - `n_terminal_length_scan.csv`: Complete results for single repeat
   - `full_length_performance.csv`: Full-length sequence performance

## Results Interpretation Guide

### Key Statistical Metrics

1. **Turning Point Identification**
   - The script identifies the first N-terminal length where performance is not significantly different from full-length
   - Criteria: Significant ratio < 50%, accuracy ≥ 90% of full-length performance

2. **Significance Ratio**
   - Proportion of repeats showing significant difference (p < 0.05) from full-length
   - High ratio (>0.7): Strong evidence of difference
   - Low ratio (<0.3): Strong evidence of similarity

3. **Performance Gap**
   - Percentage difference from full-length accuracy
   - Negative values: N-terminal performs better than full-length
   - Positive values: N-terminal performs worse

### Hypothesis Testing Regions

1. **Rapid Divergence (40-49AA)**
   - Expected: High significance ratio, large positive performance gap
   - Interpretation: Short N-terminal segments insufficient for accurate classification

2. **Transition (50-84AA)**
   - Expected: Decreasing significance ratio, decreasing performance gap
   - Interpretation: Gradual accumulation of functional information

3. **Noise Removal (85-100AA)**
   - Expected: Low significance ratio, near-zero or negative performance gap
   - Interpretation: Additional sequence removes noise without adding information

### Biological Interpretation Examples

```python
# Strong N-terminal signal:
Turning point ≤ 55AA, accuracy ratio ≥ 0.95
# Interpretation: Functional information highly concentrated in N-terminal

# Distributed signal:
Turning point ≥ 80AA, accuracy ratio < 0.90
# Interpretation: Functional information distributed across sequence
```

## Methodological Details

### ESM2 Model Configuration
- Uses `facebook/esm2_t12_35M_UR50D` (35M parameters)
- Falls back to `facebook/esm2_t6_8M_UR50D` if primary model unavailable
- Embedding extraction: Mean pooling of all residue embeddings

### SVM Configuration
- LinearSVC with automatic parameter tuning (C: 0.01, 0.1, 1.0, 10.0)
- Class weighting: Balanced for imbalanced datasets
- PCA dimensionality reduction (retains 90% variance)
- 5-fold stratified cross-validation

### Robustness Features
- Multiple SVM configurations with fallback mechanisms
- Automatic handling of rare classes (<5 samples → merged)
- PCA with randomized SVD for numerical stability
- Memory-efficient batch processing

## Troubleshooting

### Common Issues and Solutions

1. **Memory Issues**
```python
# Solution 1: Use smaller ESM2 model
# Modify line 34 in script:
model_name="facebook/esm2_t6_8M_UR50D"

# Solution 2: Reduce batch size
# Modify line 271 in script:
batch_size = 1  # Default is 2
```

2. **Slow Execution**
```python
# Use GPU acceleration
python repeat_linearSVM_Length_scan.py --device "cuda"

# Reduce number of repeats
python repeat_linearSVM_Length_scan.py --n_repeats 3
```

3. **Data Loading Errors**
```bash
# Check file formats
# FASTA: Ensure standard format with >header
# CSV: Ensure gene identifiers match FASTA headers
```

4. **Visualization Errors**
```bash
# Install required packages
pip install seaborn
# Or disable visualization in code (comment lines 1199-1369)
```

### Error Messages and Solutions

1. **"Model loading failed"**
   - Check internet connection for model download
   - Try smaller model: `facebook/esm2_t6_8M_UR50D`

2. **"No common genes found"**
   - Verify gene identifiers match between FASTA and CSV files
   - Check separator detection in CSV file

3. **"CUDA out of memory"**
   - Use CPU instead: `--device "cpu"`
   - Reduce batch size in code (line 271)

## Advanced Customization

### Modifying Analysis Parameters
```python
# In the main() function (lines 1296-1331), modify:
analyzer = RepeatedLinearSVMLengthScanAnalyzer(
    model_name="facebook/esm2_t33_650M_UR50D",  # Larger model
    device="cuda"  # Force GPU usage
)

# In robust_linear_svm_cv method (lines 342-347), adjust:
svm_configs = [
    {'C': 0.001, 'max_iter': 10000},  # More conservative
    {'C': 0.01, 'max_iter': 10000},
    {'C': 0.1, 'max_iter': 10000},
    {'C': 1.0, 'max_iter': 20000},    # More iterations
    {'C': 10.0, 'max_iter': 20000},
]
```

### Adding New Hypothesis Regions
```python
# In analyze_hypothesis_support method (lines 1109-1115), add:
hypothesis_regions = {
    'rapid_divergence': (40, 49),
    'transition': (50, 84),
    'noise_removal': (85, 100),
    'new_hypothesis': (101, 120)  # Custom region
}
```

### Extending Length Range
```python
# Command line option:
python repeat_linearSVM_Length_scan.py --min_length 20 --max_length 200

# Or modify default values in argument parser (lines 1273-1276)
parser.add_argument('--min_length', type=int, default=20)
parser.add_argument('--max_length', type=int, default=200)
```

## Citation and Acknowledgments

If using this script, please cite:

- **ESM2 Model**: Lin et al. (2022) *Nature Methods*
- **Linear SVM**: Scikit-learn library (Pedregosa et al., 2011)
- **Statistical Methods**: SciPy library (Virtanen et al., 2020)

**Methodological Reference**:  
*Repeated Linear SVM analysis for robust detection of N-terminal functional signals in plant P450 evolution*

## Support and Contact

For technical issues:
1. Check console error messages
2. Verify data file formats and locations
3. Ensure all dependencies are installed
4. Reduce analysis scale for debugging

For scientific questions about the methodology or interpretation, please consult the accompanying manuscript or contact the corresponding author.

---

**Script Version**: repeat_linearSVM_Length_scan.py v1.0  
**Last Updated**: 2025  
**Design Purpose**: Robust hypothesis testing of P450 N-terminal evolution using repeated Linear SVM analysis
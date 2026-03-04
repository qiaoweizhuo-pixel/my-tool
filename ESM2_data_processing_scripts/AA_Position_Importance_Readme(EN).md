# AA_Position_Importance.py: Enhanced Positional Importance Analysis with 6-Repeats and Positive Selection Comparison

## Script Overview
This is an enhanced version of the positional importance analysis script that implements the **"Global SVM Weight Back-Projection with Six Repeats"** method, which has been selected as the final approach for our plant P450 N-terminal signal evolution study. The script performs robust identification of key amino acid positions within the N-terminal region that are critical for subcellular localization classification.

## Method Description
The script implements the following analytical workflow:

1. **ESM2 Embeddings**: Uses the `facebook/esm2_t33_650M_UR50D` model to generate protein sequence embeddings at both global and positional levels.

2. **6-Repeat Linear SVM Analysis**:
   - Performs six independent repetitions with different random seeds (42, 123, 456, 789, 101112, 131415)
   - Each repeat includes 5-fold stratified cross-validation
   - Global SVM weights are back-projected onto individual amino acid positions
   - Positional importance scores are calculated based on absolute projection values

3. **Result Integration**:
   - Aggregates results across all six repeats
   - Calculates mean, standard deviation, median, and coefficient of variation for each position
   - Identifies consistently important positions (defined as top 20% by mean rank with CV < 0.5)

4. **Positive Selection Comparison**:
   - Compares SVM-identified key positions with known positive selection sites
   - Performs hypergeometric test for statistical significance
   - Calculates overlap ratios and Jaccard similarity coefficients

## Input Requirements

### Required Files:
1. **FASTA File**: Protein sequences in FASTA format
   - Default: `P450_unique_pep_final.fasta`
   - Should contain unique P450 protein sequences

2. **Cluster/Label File**: CSV file with gene-to-cluster mapping
   - Default: `P450_unique_pep_final_esm_clusters_genes.csv`
   - Format: Should contain at least two columns: 'Gene_List' and 'Cluster'
   - The 'Gene_List' column contains gene identifiers separated by semicolons

### Optional Configuration:
- **Positive Selection Positions**: Pre-defined list of positions under positive selection (provided as a Python list in the script)
- **N-terminal Length**: Number of amino acids to analyze from the N-terminus (default: 100)

## Dependencies

### Core Libraries:
- `torch` >= 1.10.0
- `transformers` >= 4.20.0
- `scikit-learn` >= 1.0.0
- `numpy` >= 1.21.0
- `pandas` >= 1.3.0
- `matplotlib` >= 3.5.0
- `seaborn` >= 0.11.0
- `scipy` >= 1.7.0
- `tqdm` >= 4.62.0
- `biopython` >= 1.79 (for FASTA parsing in main function)

### Installation:
```bash
pip install torch transformers scikit-learn numpy pandas matplotlib seaborn scipy tqdm biopython
```

## Usage

### Basic Execution:
```bash
python AA_Position_Importance.py
```

The script will automatically:
1. Load the FASTA and cluster files from default locations
2. Extract N-terminal sequences (first 100 amino acids)
3. Perform 6-repeat SVM analysis
4. Compare with positive selection positions
5. Generate visualizations and comprehensive reports

### Configuration Parameters (in main() function):
- `FASTA_FILE`: Path to FASTA file (default: "P450_unique_pep_final.fasta")
- `CLUSTER_FILE`: Path to cluster/label file (default: "P450_unique_pep_final_esm_clusters_genes.csv")
- `OUTPUT_DIR`: Output directory (default: "enhanced_position_importance_6repeats")
- `POSITIVE_SELECTION_POSITIONS`: List of positions under positive selection
- `N_AA`: Number of N-terminal amino acids to analyze (default: 100)
- `device`: Computation device ("cpu" or "cuda")

### Custom Analysis (Python API):
```python
from AA_Position_Importance import EnhancedPositionImportanceAnalyzer

# Initialize analyzer
analyzer = EnhancedPositionImportanceAnalyzer(model_name="facebook/esm2_t33_650M_UR50D", device="cpu")

# Run 6-repeat analysis
results = analyzer.six_repeats_analysis(sequences, labels, n_aa=100)

# Compare with positive selection positions
comparison = analyzer.compare_with_positive_selection(key_positions, positive_selection_positions)

# Generate visualizations
analyzer.visualize_integrated_results(results['integrated'], positive_selection_positions)

# Generate comprehensive report
analyzer.generate_comprehensive_report(results['integrated'], comparison)
```

## Output Structure

The script generates the following outputs:

### Main Output Directory (`enhanced_position_importance_6repeats/`):
- `comprehensive_analysis_report.txt`: Detailed analysis report with executive summary, key findings, and biological interpretations
- `integrated_position_importance.csv`: CSV file with importance scores for each position (mean, std, median, CV)
- `comparison_results.json`: JSON file with statistical comparison between SVM positions and positive selection sites
- `analysis_metadata.json`: Metadata about the analysis run

### Subdirectories:
- `plots/`: Contains all visualization files (PNG and PDF formats)
  - `integrated_position_importance.png`: Main 2x2 panel visualization
  - `svm_vs_positive_selection.png`: Comparison with positive selection positions
  - `position_stability.png`: Coefficient of variation analysis
- `repeat_*/`: Individual repeat results (6 directories, one per repeat)
  - `position_importance.csv`: Position importance scores for that repeat
  - `metadata.json`: Repeat-specific metadata (accuracy, stability metrics)

### Visualization Files:
1. **Integrated Position Importance**: Four-panel figure showing:
   - Bar plot of position importance with error bars
   - Heatmap of importance across repeats
   - Cumulative importance distribution
   - Accuracy distribution boxplot

2. **SVM vs Positive Selection**: Line plot comparing SVM importance scores with positive selection positions

3. **Position Stability**: Bar plot of coefficient of variation for each position

## Methodological Notes

### Why This Method Was Chosen (vs Local Position-Independent SVM):
1. **Global Model Consistency**: Uses a single global SVM model rather than position-specific models, ensuring consistent weight interpretation across positions
2. **Reduced Variance**: Six repeats with different random seeds provide robust estimates of position importance
3. **Statistical Rigor**: Includes comprehensive stability metrics (weight correlations, Kendall's W, coefficient of variation)
4. **Biological Relevance**: Direct comparison with positive selection data provides evolutionary context

### Key Methodological Features:
- **Weight Back-Projection**: Global SVM weights are mathematically back-projected to original amino acid positions
- **Stability Metrics**: Multiple metrics assess the reliability of position importance estimates
- **Statistical Testing**: Hypergeometric test quantifies overlap significance with positive selection
- **Comprehensive Reporting**: Generates detailed textual and visual reports

## Citation and Method Reference

When using this method in publications, please cite:

> [Your paper reference here]
>
> **Method**: "Global SVM Weight Back-Projection with Six-Repeat Analysis for Positional Importance Identification in Protein Sequences"

## Troubleshooting

### Common Issues:
1. **Memory Error**: Reduce batch size in embedding extraction or use CPU instead of GPU
2. **JSON Serialization Error**: The script includes custom serialization handlers; ensure all NumPy types are converted
3. **Missing Dependencies**: Install required libraries with `pip install -r requirements.txt`
4. **File Format Issues**: Ensure FASTA and CSV files follow the specified formats

### Performance Tips:
- Use GPU (`device="cuda"`) for faster embedding extraction
- Adjust `N_AA` parameter based on available memory
- For large datasets, consider increasing `max_length` in tokenizer or using sequence truncation

## Contact

For questions or issues related to this script, please contact:
- [Your Name/Institution]
- [Email]
- [GitHub Repository Link]

---

**Note**: This script represents the final methodological choice for the plant P450 N-terminal signal evolution study. All results in the associated manuscript should be generated using this exact implementation.
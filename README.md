# PLMA

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18857481.svg)](https://doi.org/10.5281/zenodo.18857481)
[![GitHub Release](https://github.com/qiaoweizhuo-pixel/my-tool)

**Version**: 1.0.0 | **Archived**: Zenodo DOI: [10.5281/zenodo.18857481](https://doi.org/10.5281/zenodo.18857481)  

This repository contains all the analysis code for the manuscript titled **"Protein Language Modeling and Evolutionary Analysis Reveal an N-terminal Determinant of Functional Divergence in Cytochrome P450s from Sophora. tonkinensis**.

## 📋 Repository Overview

This repository is organized to follow the **Methods section** of the manuscript in the order of execution. The workflow proceeds as:

1. **Method 2.6**: Generate evolutionary and language model features (`esm_evol_als.py`)
2. **Method 2.7.1**: Evaluate classification algorithms (`Classifier_Selection.py`)
3. **Method 2.7.2**: Optimize sequence alignment length (`repeat_linearSVM_Length_scan.py`)
4. **Method 2.8**: Identify important amino acid positions (`AA_Position_Importance.py`)

## 📁 Repository Structure

```
my-tool/
│
├──ESM_EVOL-ALS
│   ├── esm_evol_als.py           # Main entry point and workflow orchestration  Method 2.6: Feature         │   │                             # generation
│   ├── core/                     # Core computational modules
│   │   ├── __init__.py
│   │   ├── data_processor.py    # FASTA parsing and ESM2 feature extraction
│   │   ├── clustering.py        # Clustering algorithms (hierarchical + phylogenetic)
│   │   ├── metrics.py           # Evaluation metrics (V-measure, AMI)
│   │   └── utils.py             # Utility functions (distance calculations, PCoA)
│   └── visualization/           # Visualization and export modules
│       ├── __init__.py
│       ├── exporter.py          # Data export functionality (CSV, Newick)
│       └── plotter.py           # Plot generation (static and interactive)
│
├──ESM2_data_processing_scripts
│   ├── Classifier_Selection.py       # Method 2.7.1: Classifier evaluation
│   ├── Classifier_Selection.py_ReadMe(CHN).md
│   ├── Classifier_Selection.py_ReadMe(EN).md
│   ├── repeat_linearSVM_Length_scan.py  # Method 2.7.2: Length optimization
│   ├── Repeat_linearSVM_Length_scan_ReadMe(CHN).md
│   ├── Repeat_linearSVM_Length_scan_ReadMe(EN).md
│   ├── AA_Position_Importance.py     # Method 2.8: Position importance analysis
│   ├── AA_Position_Importance_Readme(CHN).md
│   └── AA_Position_Importance_Readme(EN).md
│
├── data/                             # Input Example data directory
│   ├── example.fasta                 # Protein sequences (example)
│   └── example_treefile.nwk          # Evolutionary tree corresponding to the full-length protein sequence
│
├──README(CHN)_for_ESM_EVOL_ALS.md
├──README(En)_for_ESM_EVOL_ALS.md
├──README.md
├──requirement.txt
└──LICENSE

```

## 🚀 Getting Started

### Prerequisites

- Python 3.12
- Note: Python 3.13 has not been tested with this requirement set. It may cause compatibility issues (e.g., building older packages like numpy==1.23.5 fails due to removed internal modules). We strongly recommend using Python 3.12 or lower for a smooth installation.
- Required Python packages (install via `requirements.txt`):

```bash
pip install -r requirements.txt
```

### Data Preparation

1.  Your protein sequences in FASTA format 
2. To generate an evolutionary tree based on the full-length protein sequence, it is recommended to use iqtree2 software

## 🔬 Analysis Workflow

### Step 1: Feature Generation (Method 2.6)

**Script**: `esm_evol_als.py`

This script implements **Method 2.6** of the paper. It generates the evolutionary and ESM-2 embedding features (**XXX_esm_clusters_genes.csv**) that serve as the essential input data for all downstream analyses (**Methods 2.7.1, 2.7.2, and 2.8**).

**For specific usage instructions, please refer to README(CHN)\_for_ESM_EVOL_ALS.md or README(EN)\_for_ESM_EVOL_ALS.md**

### Step 2: Classifier Evaluation (Method 2.7.1)

**Script**: `Classifier_Selection.py`

Implements **Method 2.7.1** to evaluate the performance of different machine learning classifiers for P450 protein family classification.

**For specific usage instructions, please refer to Classifier_Selection.py_ReadMe(CHN).md or Classifier_Selection.py_ReadMe(EN).md.**

### Step 3: Sequence Length Optimization (Method 2.7.2)

**Script**: `repeat_linearSVM_Length_scan.py`

Implements **Method 2.7.2** to systematically evaluate how the starting and ending of N-terminal positions of sequence affect classification performance using linear SVM.

**For specific usage instructions, please refer to Repeat_linearSVM_Length_scan_ReadMe(CHN).md or Repeat_linearSVM_Length_scan_ReadMe(EN).md.**

### Step 4: Amino Acid Position Importance (Method 2.8)

**Script**: `AA_Position_Importance.py`

Implements **Method 2.8** to identify and visualize amino acid positions that are most important for distinguishing between P450 protein families.

**For specific usage instructions, please refer to AA_Position_Importance_Readme(CHN).md or AA_Position_Importance_Readme(EN).md.**

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@software{PLMA,
  author = {Qiao Zhu},
  title = {Code for: Protein Language Modeling and Evolutionary Analysis Reveal an N-terminal Determinant of Functional Divergence in Cytochrome P450s from Sophora. tonkinensis},
  year = {2026},
  publisher = {Zenodo},
  version = {v1.0},
  doi = {DOI 10.5281/zenodo.18857481.},
  url = {https://doi.org/10.5281/zenodo.18857481}
}
```

## 🤝 Contributing

For questions or issues regarding the code, please open an issue on GitHub.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

**Note**: This code is associated with the manuscript currently under review. Please contact the corresponding author for pre-print access.

## 📜 Version History

- **v1.0.0** (2024-03-04): Initial release corresponding to manuscript submission
  - Archived at Zenodo: [10.5281/zenodo.18857481](https://doi.org/10.5281/zenodo.18857481)
  - Includes all four core analysis scripts

**Zenodo (archived):** v1.0.0 – DOI: 10.5281/zenodo.18857481
This is the fixed version cited in the paper, ensuring exact reproducibility.

**GitHub (development):** May contain bug fixes and updates
For the latest features, please refer to the most recent commits, but note that there might be minor discrepancies compared to the results in the paper.

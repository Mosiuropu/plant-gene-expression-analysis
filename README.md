# 🧬 Plant Gene Expression Analysis Under Abiotic Stress

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

A comprehensive bioinformatics pipeline for analyzing **Arabidopsis thaliana** gene expression data under abiotic stress conditions (Drought, Heat, and Cold). This project demonstrates reproducible research workflows in plant genomics using Python.

---

## 📋 Table of Contents

- [Project Overview](#-project-overview)
- [Key Features](#-key-features)
- [Repository Structure](#-repository-structure)
- [Installation](#-installation)
- [Usage](#-usage)
- [Methodology](#-methodology)
- [Results](#-results)
- [Biological Significance](#-biological-significance)
- [Technical Skills Demonstrated](#-technical-skills-demonstrated)
- [License](#-license)

---

## 🎯 Project Overview

Plants are sessile organisms that must constantly adapt to environmental stresses. Understanding the **transcriptional reprogramming** that occurs under stress conditions is critical for developing climate-resilient crops.

This project simulates *Arabidopsis thaliana* gene expression data under four conditions:
- **Control** — Optimal growth conditions
- **Drought** — Water deficit stress
- **Heat** — Elevated temperature stress
- **Cold** — Low temperature stress

The pipeline performs end-to-end analysis from raw expression data to biological interpretation, including quality control, statistical analysis, and publication-quality visualizations.

---

## ✨ Key Features

- **Realistic Data Simulation** — Generates expression data with known stress-marker genes and biologically meaningful patterns
- **Quality Control** — Filters lowly-expressed genes, detects outlier samples, and normalizes data
- **Dimensionality Reduction** — PCA with publication-quality scatter plots
- **Differential Expression** — Welch's t-test with Benjamini-Hochberg FDR correction
- **GO Enrichment** — Fisher's exact test for biological process overrepresentation
- **Publication-Quality Figures** — PCA plots, volcano plots, heatmaps, boxplots, GO bar charts
- **Reproducible** — Fully scripted pipeline with random seed control
- **Interactive Notebook** — Jupyter notebook for exploratory data analysis

---

## 📁 Repository Structure

```
plant-gene-expression-analysis/
│
├── src/                          # Core Python modules
│   ├── __init__.py               # Package initializer
│   ├── data_generation.py        # Simulated expression data generation
│   ├── preprocessing.py          # QC filtering, normalization, outlier detection
│   ├── analysis.py               # PCA, differential expression, GO enrichment
│   └── visualization.py          # All plotting functions
│
├── notebooks/
│   └── 01_exploratory_analysis.ipynb  # Interactive Jupyter notebook
│
├── scripts/
│   └── run_analysis.py           # Main pipeline orchestration script
│
├── data/                         # Generated expression data (CSV)
│   ├── expression_data.csv       # Gene expression matrix (genes × samples)
│   └── metadata.csv             # Sample metadata with condition labels
│
├── results/                      # Pipeline output
│   └── figures/                  # Generated figures (PNG)
│       ├── 01_pca_plot.png
│       ├── 02_volcano_*.png
│       ├── 03_heatmap_*.png
│       ├── 04_marker_expression.png
│       └── 05_go_enrichment_*.png
│
├── requirements.txt              # Python dependencies
├── LICENSE                       # MIT License
└── README.md                     # This file
```

---

## 🔧 Installation

### Prerequisites

- **Python 3.9+** installed on your system
- Basic familiarity with the terminal/command line

### Setup Instructions

```bash
# 1. Clone the repository
git clone https://github.com/Mosiuropu/plant-gene-expression-analysis.git
cd plant-gene-expression-analysis

# 2. Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate        # On Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Verify installation
python -c "import numpy; import pandas; print('All dependencies installed!')"
```

---

## 🚀 Usage

### Quick Start (Full Pipeline)

Run the complete analysis pipeline with a single command:

```bash
python scripts/run_analysis.py
```

This will:
1. Generate simulated expression data for 5,000 Arabidopsis genes
2. Filter and normalize the data
3. Run PCA, differential expression, and GO enrichment
4. Generate all figures
5. Create a markdown analysis report

### Command-Line Options

```bash
python scripts/run_analysis.py --help

# Custom data directory
python scripts/run_analysis.py --data-dir my_data --output-dir my_results

# Use existing data (skip generation)
python scripts/run_analysis.py --data-dir data --no-generate

# Different random seed
python scripts/run_analysis.py --seed 123
```

### Interactive Analysis (Jupyter)

```bash
jupyter notebook notebooks/01_exploratory_analysis.ipynb
```

### Module Usage (Python API)

```python
from src.data_generation import generate_dataset
from src.preprocessing import run_preprocessing
from src.analysis import run_pca, differential_expression, go_enrichment

# Generate data
expression_df, metadata_df = generate_dataset()

# Preprocess
preprocessed = run_preprocessing(expression_df, metadata_df)

# Run PCA
pca_results = run_pca(preprocessed['expression'])

# Differential expression
de_results = differential_expression(
    preprocessed['expression'], 'Control', 'Drought',
    preprocessed['metadata']
)

# GO enrichment
go_results = go_enrichment(de_results, expression_df.index.tolist())
```

---

## 🔬 Methodology

### Data Simulation

The expression data is designed to mimic real RNA-seq data from *Arabidopsis thaliana*:

| Component | Details |
|-----------|---------|
| **Organism** | *Arabidopsis thaliana* (model plant) |
| **Genes** | 5,000 (subset for computational efficiency) |
| **Conditions** | Control, Drought, Heat, Cold |
| **Replicates** | 4 per condition (16 total samples) |
| **Expression units** | log₂(TPM+1) |
| **Marker genes** | 25 well-characterized stress-responsive genes |

### Statistical Methods

| Analysis | Method | Reference |
|----------|--------|-----------|
| **Normalization** | Quantile normalization | Bolstad et al., 2003 |
| **Dimensionality reduction** | Principal Component Analysis | Pearson, 1901 |
| **Differential expression** | Welch's t-test | Welch, 1947 |
| **Multiple testing correction** | Benjamini-Hochberg FDR | Benjamini & Hochberg, 1995 |
| **GO enrichment** | Fisher's exact test | Fisher, 1922 |

### Known Stress Marker Genes

| Gene | Function | Drought | Heat | Cold |
|------|----------|:-------:|:----:|:----:|
| `RD29A` | Desiccation protection | ↑↑↑ | ↑ | ↑↑↑ |
| `HSP70` | Protein folding chaperone | ↑ | ↑↑↑ | — |
| `CBF3` | Cold acclimation master regulator | — | ↓ | ↑↑↑ |
| `NCED3` | ABA biosynthesis (drought signaling) | ↑↑↑ | — | ↑ |
| `COR15A` | Membrane stabilization (cold) | — | — | ↑↑↑ |
| `RBCS1A` | Photosynthesis (Rubisco) | ↓↓ | ↓↓ | ↓ |
| `WRKY33` | Multi-stress transcription factor | ↑↑ | ↑↑ | ↑ |

---

## 📊 Results

After running the pipeline, you will find in `results/figures/`:

| Figure | Description |
|--------|-------------|
| `01_pca_plot.png` | PCA separation of samples by condition |
| `02_volcano_Drought_vs_Control.png` | DEGs under drought stress |
| `02_volcano_Heat_vs_Control.png` | DEGs under heat stress |
| `02_volcano_Cold_vs_Control.png` | DEGs under cold stress |
| `03_heatmap_Drought_vs_Control.png` | Expression heatmap of top DEGs |
| `04_marker_expression.png` | Known stress marker gene expression |
| `05_go_enrichment_Drought_vs_Control.png` | GO enrichment — Drought vs Control |
| `05_go_enrichment_Heat_vs_Control.png` | GO enrichment — Heat vs Control |
| `05_go_enrichment_Cold_vs_Control.png` | GO enrichment — Cold vs Control |

A complete **analysis report** is generated at `results/analysis_report.md`.

---

## 🌱 Biological Significance

Understanding plant stress responses has critical implications:

- **Food Security**: Climate change threatens global crop yields. Understanding stress response mechanisms enables breeding of stress-tolerant varieties.
- **Conservation**: Predicting how natural plant populations will respond to changing climates.
- **Biotechnology**: Identifying target genes for genetic engineering of stress-resistant crops.

The CBF (C-Repeat Binding Factor) regulon, ABA (Abscisic Acid) signaling pathway, and heat shock protein networks are among the most well-studied stress response systems in plants — all of which are represented in this analysis.

---

## 💻 Technical Skills Demonstrated

| Skill | Evidence |
|-------|----------|
| **Python Programming** | Modular OOP design, type hints, docstrings |
| **Data Analysis** | NumPy, pandas, SciPy, scikit-learn |
| **Statistics** | Hypothesis testing, multiple testing correction |
| **Bioinformatics** | Gene expression analysis, GO enrichment |
| **Data Visualization** | matplotlib, seaborn |
| **Git/GitHub** | Version control, meaningful commits, README |
| **Reproducibility** | Seed-based simulation, automated pipeline |
| **Jupyter** | Interactive notebooks for exploratory analysis |
| **Documentation** | Comprehensive README, docstrings, comments |
| **Project Structure** | Modular package layout, separation of concerns |

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  <sub>Created for research demonstration purposes | Arabidopsis stress transcriptomics</sub>
</p>

# ADNeuronNet

This repository contains the code for data analysis and figure generation for the manuscript:

**"Interactome mapping in human excitatory neurons reveals novel drivers in Alzheimer's disease."**

## Project Description

ADNeuronNet provides reproducible scripts for the statistical analyses and visualizations presented in our study of protein–protein interaction (interactome) networks in human excitatory neurons, with a focus on identifying novel molecular drivers of Alzheimer's disease. The analyses span multi-sample proportion comparisons, volcano-plot visualization of differential interactors, and pairwise statistical testing of protein abundance between wild-type and mutant conditions.

## Repository Structure

```
ADNeuronNet/
├── Fig2/
│   ├── Fig2.ipynb          # Jupyter notebook for Figures 2B & 2C
│   └── Fig2_data.xlsx      # Source data (neuron-expressed & neuron-specific interactors)
├── Fig4/
│   ├── Fig4a.R             # R script for the BIN1 network volcano plot
│   └── Fig4_data.xlsx      # Source data (BIN1v1 network)
├── Fig5/
│   └── Fig5g.ipynb         # Jupyter notebook for RIN3 WT vs. W63C comparison
├── LICENSE
└── README.md
```

## Figures

| Figure | Script | Description |
|--------|--------|-------------|
| **Fig 2B–C** | `Fig2/Fig2.ipynb` | Bar plots comparing the percentage of neuron-expressed and neuron-specific interactors across conditions, with pairwise binomial proportion tests and significance brackets. |
| **Fig 4A** | `Fig4/Fig4a.R` | Enhanced volcano plot of the BIN1 interactome network, highlighting protein complexes (APC/C, RIN, AP2, CLT, PIK3, TIM) with color-coded labels. |
| **Fig 5G** | `Fig5/Fig5g.ipynb` | Bar plot with individual data points comparing normalized intensity of RIN3 wild-type vs. RIN3 W63C mutant, using an independent-samples t-test. |

## Requirements

### Python (Fig 2, Fig 5)

- Python ≥ 3.8
- numpy
- pandas
- scipy
- seaborn
- matplotlib

Install with pip:

```bash
pip install numpy pandas scipy seaborn matplotlib openpyxl
```

### R (Fig 4)

- R ≥ 4.0
- tidyverse
- readxl
- EnhancedVolcano
- ggrepel
- RColorBrewer
- dplyr

Install in R:

```r
install.packages(c("tidyverse", "readxl", "ggrepel", "RColorBrewer"))
# EnhancedVolcano is available from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
```

## Usage

### Jupyter Notebooks (Fig 2 & Fig 5)

Open and run all cells in order:

```bash
cd Fig2 && jupyter notebook Fig2.ipynb
cd Fig5 && jupyter notebook Fig5g.ipynb
```

### R Script (Fig 4)

Run the script from the `Fig4/` directory so that relative data paths resolve correctly:

```bash
cd Fig4 && Rscript Fig4a.R
```

## License

This project is licensed under the [MIT License](LICENSE).

Copyright © 2026 Yu Sun
# ADNeuronNet

This repository contains the code for data analysis and figure generation for the manuscript:

**"Interactome mapping in human excitatory neurons reveals novel drivers in Alzheimer's disease."**

Developed by [Yu Sun](https://github.com/YuSun795) and [Tianyu Xia](https://github.com/TianyuXia123) form Haiyuan Yu's Lab.

## Project Abstract

Alzheimer’s disease (AD) is an irreversible neurodegenerative disease defined by its molecular hallmarks - amyloid beta peptide plaques and neurofibrillary Tau tangles. Despite significant progress that has been made in uncovering a large number of genetic risk factors through extensive genomic sequencing and genetic studies, the molecular mechanisms driving AD-associated pathology and cognitive decline remain poorly understood. Therefore, alongside the identification of more risk genes, it is also paramount to study how these genes function and influence each other within the cellular pathways and overall molecular networks in AD-relevant brain cell types. However, current human protein-protein interactome datasets were all generated in either yeast or generic human cell lines. Consequently, many important neuronal interactions, especially neuron-specific ones, have yet been discovered. To address this critical gap, we developed a highly scalable, high-quality interactome mapping pipeline in human excitatory neurons derived from induced pluripotent stem cells (iPSC), and generated a comprehensive, neuron-specific interactome map, named ADNeuronNet, for key AD risk genes. ADNeuronNet consists of 1,767 high-confidence interactions among 1,189 proteins and is the only dataset enriched with neuron-specific genes when compared to known protein interactions, including previous large-scale interactome maps, for the same baits in the literature. Within ADNeuronNet, we identified 1,375 novel interactions, many of which are likely neuron specific. For example, we identified a neuron-specific interactor, RIN2, for major AD risk factor BIN1 and confirmed RIN2’s function in recruiting BIN1 to early endosomes, a process that has been well-associated with AD etiology. Additionally, we performed quantitative interaction perturbation analyses on AD risk genes with AD-associated mutations or isoforms and identified significant changes in 91 protein interactions among 11 different protein variants. Finally, we found that subunits from the anaphase-promoting complex/cyclosome (APC/C), other novel BIN1 interactors identified by ADNeuronNet, mediated modulation of Tau-aggregation in neurons via regulation of APOE expression, uncovering a previously unrecognized BIN1-APC/C-APOE regulatory axis in AD pathobiology. Overall, these findings illustrate how our neuron-specific ADNeuronNet can be leveraged to uncover new risk gene candidates and cellular pathways that help advance our understanding of molecula

## Content in this repository

### `data_processing/`

- **`process_af3_results.py`** — Processes AlphaFold3 (AF3) protein structure prediction results to extract quality metrics (ipTM, pTM, ranking score) and map AD-associated missense mutations onto predicted protein-protein interaction interfaces. Outputs summary score tables and mutation-to-interface mapping files.

- **`process_literature_interactomes.ipynb`** — Consolidates and standardizes three major literature-derived PPI databases (HINT, BioPlex, OpenCell) with HGNC-approved gene symbols. Maps UniProt IDs to approved gene symbols and produces a unified literature PPI reference set.

### `Fig2/`

- **`Fig2_data.xlsx`** — Source data workbook containing multiple sheets used to generate Figures 2a–h.

- **`Fig2a-e.ipynb`** — Analyzes PPI enrichment in ADNeuronNet against literature databases, characterizes neuron-specific interactors, and evaluates CRISPRi hit overlap. Generates Figures 2a–e using Fisher's exact tests with multiple testing correction.

- **`Fig2f_enrichr_pathway_all_Jan2026.py`** — Performs pathway enrichment analysis on ADNeuronNet interactors using hypergeometric testing against KEGG/GO pathway databases, with FDR-BH correction.

- **`Fig2f_pwyenrich_plot.R`** — Visualizes pathway enrichment results as a publication-quality bubble plot (top 10 significant pathways by adjusted p-value, colored by enrichment ratio).

- **`Fig2g_enrich_HPA.R`** — Creates a protein/tissue enrichment bubble plot from Human Protein Atlas (HPA) GOSt database annotations for ADNeuronNet interactors.

- **`Fig2h_GSEA_permutation_analysis.py`** — Performs Gene Set Enrichment Analysis (GSEA) with degree-preserving permutation testing to assess ADNeuronNet interactor enrichment across four endophenotypes (Tau, Lysosome, Synapse, Mitochondrial).

- **`Fig2h_Permutation_Endophenotype_nature_histogram_violin.R`** — Creates a nature-style violin plot comparing observed vs. permuted odds ratios across the four endophenotypes, with empirical p-value annotations.

### `Fig3/`

- **`Fig3e.ipynb`** — Validates ADNeuronNet PPI quality by comparing AlphaFold3 ipTM score distributions across positive controls, negative controls, random pairs, and literature interactions using Mann-Whitney U tests.

### `Fig4/`

- **`Fig4a.R`** — Generates an enhanced volcano plot for BIN1 network differential expression analysis, highlighting APC/C complex members and other key interactors using the EnhancedVolcano package.

- **`Fig4c.ipynb`** — Displays RNA expression levels (TPM) of RIN2 and RIN3 across human cell lines (HCT116, HEK293) and brain regions (amygdala, cerebellum, cerebral cortex, hippocampal formation, thalamus) as an annotated heatmap using HPA data.

### `Fig5/`

- **`Fig5g.ipynb`** — Performs a statistical comparison (two-sample t-test) of normalized protein intensity between RIN3 wild-type and W63C mutant, visualized as a bar plot with individual data points and significance annotation.

### `Fig6/`

- **`ANAPC realted R script .R`** — Comprehensive analysis of ANAPC2 knockdown vs. scrambled control: generates a volcano plot of differentially expressed genes, a GO enrichment bar plot for upregulated genes, and a GO enrichment dot plot for downregulated genes.

- **`ANAPC2_shc_Free_filtered_results.csv`** — Differential gene expression results (log2 fold changes, adjusted p-values) from ANAPC2 knockdown experiments. Serves as input data for the ANAPC-related R script.

## License

This project is licensed under the [MIT License](LICENSE).
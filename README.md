# Genomic basis of granulovirus resistance in *Plodia interpunctella*

This repository contains data and analysis code associated with the manuscript:

**"Genomic structure, not life-history strategy, predicts granulovirus resistance in *Plodia interpunctella*"**

Signe White & Mike Boots, University of California Berkeley

---

## Repository structure

### `R_code/`
R scripts for all analyses reported in the manuscript, in order of execution:

| Script | Description |
|--------|-------------|
| `01_life_history.R` | Life-history PCA, LD50 calculation, dose-response GLMs, Figures 1–2 and Supplementary Figures S1–S6 |
| `02_genomic_pca.R` | Genomic PCA, neighbour-joining tree, PC1 vs LD50 correlation, Figures 3–4 and Supplementary Figures S7–S8 |
| `03_manhattan_wza.R` | Sliding-window PCA-LD50 correlation, WZA, Manhattan plots, QQ plots, Figure 5 and Supplementary Figure S9 |
| `04_clustering.R` | K-means clustering of LD50 values, Supplementary Figures S10–S12 |
| `05_ROH_figure.R` | Runs of homozygosity analysis and figure, Supplementary Figure S13 |

Scripts should be run in order. `01_life_history.R` must be run first as it generates `ld50_calculated.csv`, which is used by all subsequent scripts.

### `data/`
Input data files required to run the analysis scripts:

| File | Description |
|------|-------------|
| `Pupal_Weights_Days_to_Development_Growth_Rates.csv` | Individual-level life-history measurements for all populations |
| `Infection_Growth.csv` | Individual-level infection assay data |
| `Binomial_Mortality_Data.csv` | Binary infection outcomes by dose and population, used for LD50 calculation |
| `ld50_calculated.csv` | Population-mean LD50 values calculated by `01_life_history.R` |
| `stock_pca_results.csv` | Genomic PCA scores for inbred lines and stock population |
| `genetic_distances.csv` | Pairwise genetic distance matrix among inbred lines |
| `snprelate_1mb.csv` | Per-window genomic PC1 scores for all inbred lines, used for Manhattan plot and WZA analyses |

---

## Software requirements

All analyses were run in R. Key package versions:

- `MASS` — LD50 estimation via `dose.p`
- `vcfR` v1.16.0 — VCF file handling
- `adegenet` v2.1.11 — genomic PCA
- `ape` v5.8 — neighbour-joining tree
- `SNPRelate` v1.38 — sliding-window PCA
- `ggplot2`, `patchwork`, `ggrepel` — figures

Computationally intensive analyses (sliding-window PCA, WZA) were run on the Savio HPC cluster at UC Berkeley (`03_manhattan_wza.R`).

---

## Raw genomic data

Raw sequencing reads (FASTQ files) for all 12 inbred lines and the stock population will be deposited in the NCBI Sequence Read Archive upon manuscript acceptance.

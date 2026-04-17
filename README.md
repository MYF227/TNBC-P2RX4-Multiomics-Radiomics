# P2RX4 in Triple-Negative Breast Cancer — Analysis Code

This repository contains the complete R analysis code accompanying the manuscript:

> **"Preliminary Integrative Analysis Reveals P2RX4-Driven Immune Remodeling and Peritumoral Phenotypes in TNBC"**  
> *(Authors, Journal, Year — update before submission)*

---

## Repository Structure

```
.
├── config.R                        # ★ START HERE — set all data paths
├── Figure2_scRNAseq_Analysis.R     # Single-cell UMAP + P2RX4 feature plots
├── Figure3_Spatial_M2_Correlation.R# Spatial transcriptomics vs CD163
├── Figure4_TCGA_Correlation.R      # TCGA-BRCA bulk RNA correlation panel
├── Figure5_GSEA_Hypoxia.R          # GSEA Hallmark Hypoxia enrichment
├── Figure6_Subtype_Specificity.R   # Pan-BRCA vs TNBC subtype comparison
├── Figure7_KM_Survival.R           # Kaplan-Meier OS analysis (Pan-BRCA)
├── Figure8A_Radiomics_Correlation.R# P2X4 vs radiomics feature scatter
├── Figure9_Radiomics_Heatmap.R     # Peritumoral habitat ratio heatmap
├── output/                         # Auto-created; all figures saved here
└── README.md
```

---

## Data Availability
Data Availability
The raw data are publicly available via GEO (GSE176078), TCGA, and TCIA.
Note on Processed Data: To protect the integrity of the peer-review process, the intermediate processed datasets (e.g., Seurat objects, extracted radiomics feature matrices) are currently not publicly shared in this repository. However, all intermediate data required to fully reproduce the results are available to reviewers upon reasonable request to the corresponding author during the review process. Full datasets will be made strictly open-access upon the formal acceptance and publication of this manuscript.
| Dataset | Source | Accession / URL |
|---|---|---|
| Single-cell RNA-seq (TNBC & ER+) | GEO | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) |
| TCGA-BRCA bulk RNA-seq (TPM) | UCSC Xena | [TCGA-BRCA.star_tpm](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.star_tpm&host=https%3A%2F%2Fgdc.xenahubs.net) |
| TCGA-BRCA clinical data | UCSC Xena | [TCGA-BRCA.clinical.tsv](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net) |
| Radiomics features | TCIA | [TCIA BRCA Collection](https://www.cancerimagingarchive.net/) |

> **Note**: Raw data files are **not** included in this repository due to size and access restrictions. Please download from the sources above and update the paths in `config.R`.

---

## Quick Start

### Step 1 — Install required R packages

Run the following once in R/RStudio:

```r
# CRAN packages
install.packages(c(
  "Seurat", "Matrix", "ggplot2", "patchwork", "ggpubr",
  "data.table", "dplyr", "tidyr", "tidyverse",
  "survival", "survminer", "pheatmap", "RColorBrewer"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "enrichplot", "msigdbr"))
```

### Step 2 — Configure paths

Open `config.R` and set each variable to match your local directory layout:

```r
SCRNA_DATA_DIR  <- "path/to/GSE176078_RAW"
TCGA_DATA_DIR   <- "path/to/TCGA_BRCA"
TPM_FILE_PATH   <- "path/to/TCGA-BRCA.star_tpm"
TNBC_COUNTS_CSV <- "path/to/TNBC_66_counts_matrix.csv"
CLINICAL_TSV    <- "path/to/TCGA-BRCA.clinical.tsv"
RADIOMICS_CSV   <- "path/to/P2X4_Radiomics_LongFormat.csv"
OUTPUT_DIR      <- "output"
```

### Step 3 — Run figure scripts

Each script is self-contained. Source `config.R` automatically and save output to `OUTPUT_DIR`:

```r
source("Figure2_scRNAseq_Analysis.R")
source("Figure4_TCGA_Correlation.R")
# ... etc.
```

---

## Figure-by-Figure Notes

### Figure 2 — scRNA-seq UMAP & P2RX4 expression
- Processes 10x Genomics sparse matrices from GSE176078 samples CID44991 (TNBC) and CID4535 (ER+).  
- Runs standard Seurat workflow: normalisation → variable features → PCA → UMAP → clustering.  
- Cell types are annotated using a curated `name_map` based on marker gene DotPlots.

### Figure 3 — Spatial transcriptomics correlation
- Requires a pre-built Seurat spatial object (`st_obj`) from Visium data for sample CID44971.  
- See the **PREREQUISITE** block at the top of `Figure3_Spatial_M2_Correlation.R` for loading instructions.

### Figure 4 — TCGA bulk RNA correlations
- Matches TNBC barcodes by truncating 15-digit TCGA IDs to 12 digits.  
- Gene lookup uses both Ensembl ID (primary) and HGNC symbol (fallback).

### Figure 5 — GSEA Hypoxia enrichment
- Requires upstream objects `counts_clean`, `group`, and `m_t2g_final`.  
- `m_t2g_final` can be obtained via `msigdbr` (see script header for one-liner).

### Figure 6 — Subtype specificity (P2RX4 vs VIM)
- Demonstrates that the negative P2RX4–VIM correlation is TNBC-specific.  
- Random seed fixed at `set.seed(2026)` for reproducibility of control sampling.

### Figure 7 — Kaplan-Meier survival (Pan-BRCA)
- Uses optimal cutpoint (`surv_cutpoint`) for P2RX4 grouping.  
- Exports TIFF (journal submission) + EPS (vector) + PDF (preview).

### Figure 8A — Radiomics scatter
- 4 pilot cases only; requires `GTV_df` radiomics wide-format table.  
- See **PREREQUISITE** block in script for guidance.

### Figure 9 — Radiomics heatmap
- Computes Habitat_3mm / GTV ratio for each radiomics feature.  
- Top features selected by Spearman correlation with P2X4 TPM.  
- `pivot_wider()` is used (modern tidyr replacement for deprecated `spread()`).

---

## Reproducibility

Each script prints `sessionInfo()` at completion. The key package versions used to generate the published figures are listed below:

| Package | Version |
|---|---|
| R | 4.5.2 |
| Seurat | 5.x |
| ggplot2 | 3.5.x |
| clusterProfiler | 4.x |
| survminer | 0.4.9 |
| data.table | 1.15.x |
| patchwork | 1.2.x |

> Run `sessionInfo()` in your R session and compare if results differ.

---

## Citation

If you use this code, please cite the associated manuscript (citation TBD upon acceptance).

---

## License

This code is released under the **MIT License**. See `LICENSE` for details.

---

## Contact

For questions about the analysis code, please open a GitHub Issue or contact the corresponding author.

# =============================================================================
# config.R — Unified Path & Parameter Configuration
# Project: P2RX4 in Triple-Negative Breast Cancer
#
# ► INSTRUCTIONS FOR USERS:
#   1. Edit ONLY this file before running any figure script.
#   2. Set each path to match your local data directory.
#   3. All figure scripts will source() this file automatically.
# =============================================================================

# --------------------------------------------------------------------------
# 1. Root directories (modify these to match your environment)
# --------------------------------------------------------------------------

# Where the 10x Genomics single-cell raw data folders are stored
SCRNA_DATA_DIR  <- "data/scRNA"          # contains GSM5354532_CID44991/ etc.

# Where the TCGA-BRCA bulk RNA-seq files are stored
TCGA_DATA_DIR   <- "data/TCGA_BRCA"

# Where the TCGA-BRCA TPM matrix file is located (full path to file)
TPM_FILE_PATH   <- "data/TCGA_BRCA/TCGA-BRCA.star_tpm"   # may end in .tsv or .txt

# Where the TNBC 66-sample count matrix CSV is
TNBC_COUNTS_CSV <- "data/TCGA_BRCA/TNBC_66_counts_matrix.csv"

# Where TCGA clinical data is stored
CLINICAL_TSV    <- "data/TCGA_BRCA/TCGA-BRCA.clinical.tsv"

# Where radiomics data (long-format CSV) is stored
RADIOMICS_CSV   <- "data/Radiomics/P2X4_Radiomics_LongFormat.csv"

# --------------------------------------------------------------------------
# 2. Output directory (all figures will be saved here)
# --------------------------------------------------------------------------

OUTPUT_DIR <- "output"

# --------------------------------------------------------------------------
# 3. Auto-create output directory if it does not exist
# --------------------------------------------------------------------------

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

message("✅ config.R loaded. Output will be saved to: ", normalizePath(OUTPUT_DIR, mustWork = FALSE))

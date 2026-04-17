# =============================================================================
# Figure 4 — Expression correlation between P2RX4 and key molecules in TNBC
#
# Data source : TCGA-BRCA TPM matrix (UCSC Xena / GDC portal)
#               TNBC sample list derived from TNBC_66_counts_matrix.csv
# Output      : Figure4_TNBC_Correlations.tiff / .eps
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ---------- 1. Load TNBC sample list (12-digit TCGA barcodes) -----------------
counts_data  <- read.csv(TNBC_COUNTS_CSV, stringsAsFactors = FALSE)
tnbc_samples <- gsub("\\.", "-", colnames(counts_data)[-c(1, 2)])

# ---------- 2. Load TPM matrix ------------------------------------------------
message("Loading TPM matrix (this may take a moment)...")

# Auto-detect file extension
tpm_file <- TPM_FILE_PATH
if (!file.exists(tpm_file)) tpm_file <- paste0(tpm_file, ".tsv")
if (!file.exists(tpm_file)) tpm_file <- sub("\\.tsv$", ".txt", tpm_file)
if (!file.exists(tpm_file)) stop("TPM file not found. Check TPM_FILE_PATH in config.R")

tpm_raw    <- fread(tpm_file, data.table = FALSE)
rownames(tpm_raw) <- tpm_raw[, 1]
tpm_matrix <- tpm_raw[, -1]
rm(tpm_raw); gc()

# Standardise column names to 12-digit barcodes
colnames(tpm_matrix) <- substr(gsub("\\.", "-", colnames(tpm_matrix)), 1, 12)

# Intersect with TNBC sample list
valid_samples <- intersect(tnbc_samples, colnames(tpm_matrix))
message("✅ TNBC samples matched: ", length(valid_samples))
if (length(valid_samples) == 0) stop("No samples matched. Check barcode format in config.R files.")

# ---------- 3. Gene extraction helper -----------------------------------------
# NOTE: Ensembl IDs below match the TCGA-BRCA star_tpm file used in this study.
#       Verify against your own TPM file header if IDs differ.
get_gene_exp <- function(symbol, ensg, pdl1_alias = NULL) {
  idx <- grep(ensg,    rownames(tpm_matrix), value = TRUE)[1]
  if (is.na(idx)) idx <- grep(paste0("^", symbol, "(\\.|$|\\|)"),
                               rownames(tpm_matrix), value = TRUE)[1]
  if (is.na(idx) && !is.null(pdl1_alias))
    idx <- grep(paste0("^", pdl1_alias, "(\\.|$|\\|)"),
                rownames(tpm_matrix), value = TRUE)[1]
  if (is.na(idx)) {
    message("WARNING: Gene not found: ", symbol, " — returning NA")
    return(rep(NA_real_, length(valid_samples)))
  }
  as.numeric(unlist(tpm_matrix[idx, valid_samples]))
}

# ---------- 4. Build expression data frame ------------------------------------
df <- data.frame(
  P2RX4 = get_gene_exp("P2RX4", "ENSG00000258982"),
  P2RX7 = get_gene_exp("P2RX7", "ENSG00000089041"),
  PD_L1 = get_gene_exp("PD-L1", "ENSG00000120217", pdl1_alias = "CD274"),
  VIM   = get_gene_exp("VIM",   "ENSG00000026025"),
  CDH1  = get_gene_exp("CDH1",  "ENSG00000039068"),
  VEGFA = get_gene_exp("VEGFA", "ENSG00000112715"),
  HIF1A = get_gene_exp("HIF1A", "ENSG00000100644")
)
message("Data frame head:"); print(head(df))

# ---------- 5. Plotting function ----------------------------------------------
plot_cor <- function(data, x_col, y_col, pt_color, plot_title) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(color = pt_color, alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black",
                fill = pt_color, alpha = 0.2, linewidth = 1) +
    stat_cor(method = "spearman",
             label.x.npc = 0.05, label.y.npc = 0.95, size = 5) +
    labs(title = plot_title,
         x = paste(x_col, "Expression (TPM)"),
         y = paste(y_col, "Expression (TPM)")) +
    theme_classic(base_size = 14) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title  = element_text(face = "bold", size = 14),
      axis.text   = element_text(color = "black", size = 12),
      axis.line   = element_line(linewidth = 0.8, color = "black"),
      axis.ticks  = element_line(linewidth = 0.8, color = "black"),
      plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
    )
}

pA <- plot_cor(df, "P2RX4", "PD_L1", "#808080", "P2RX4 vs PD-L1 in TNBC")
pB <- plot_cor(df, "P2RX7", "PD_L1", "#984EA3", "P2RX7 vs PD-L1 in TNBC")
pC <- plot_cor(df, "P2RX4", "VIM",   "#E41A1C", "P2RX4 vs VIM in TNBC")
pD <- plot_cor(df, "P2RX4", "CDH1",  "#4DAF4A", "P2RX4 vs CDH1 in TNBC")
pE <- plot_cor(df, "P2RX4", "VEGFA", "#A65628", "P2RX4 vs VEGFA in TNBC")
pF <- plot_cor(df, "P2RX4", "HIF1A", "#377EB8", "P2RX4 vs HIF1A in TNBC")

final_plot <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

# ---------- 6. Save -----------------------------------------------------------
file_base <- file.path(OUTPUT_DIR, "Figure4_TNBC_Correlations")

ggsave(paste0(file_base, ".tiff"), plot = final_plot,
       width = 12, height = 15, dpi = 300, compression = "lzw")
ggsave(paste0(file_base, ".eps"),  plot = final_plot,
       device = "eps", width = 12, height = 15)

message("✅ Figure 4 saved to: ", OUTPUT_DIR)

sessionInfo()

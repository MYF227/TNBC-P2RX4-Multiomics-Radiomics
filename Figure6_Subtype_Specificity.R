# =============================================================================
# Figure 6 — Subtype-specific analysis of P2RX4 and VIM correlation
#
# Data source : TCGA-BRCA TPM matrix; TNBC 66-sample list
# Output      : Figure6_Subtype_Specificity.tiff / .pdf
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(data.table)
library(tidyr)       # pivot_wider (replaces deprecated spread())
library(dplyr)
library(ggplot2)
library(ggpubr)

# ---------- 1. Load TPM matrix ------------------------------------------------
message("Loading TPM matrix...")

tpm_file <- TPM_FILE_PATH
if (!file.exists(tpm_file)) tpm_file <- paste0(tpm_file, ".tsv")
if (!file.exists(tpm_file)) tpm_file <- sub("\\.tsv$", ".txt", tpm_file)
if (!file.exists(tpm_file)) stop("TPM file not found. Check TPM_FILE_PATH in config.R")

tpm_raw    <- fread(tpm_file, data.table = FALSE)
rownames(tpm_raw) <- tpm_raw[, 1]
tpm_matrix <- tpm_raw[, -1]
rm(tpm_raw); gc()

colnames(tpm_matrix) <- substr(gsub("\\.", "-", colnames(tpm_matrix)), 1, 12)

# ---------- 2. Define sample groups -------------------------------------------
counts_data    <- read.csv(TNBC_COUNTS_CSV, stringsAsFactors = FALSE)
tnbc_samples   <- substr(gsub("\\.", "-", colnames(counts_data)[-c(1, 2)]), 1, 12)
valid_samples  <- intersect(tnbc_samples, colnames(tpm_matrix))

all_unique_samples <- unique(colnames(tpm_matrix))
non_tnbc_samples   <- setdiff(all_unique_samples, valid_samples)

set.seed(2026)
control_samples <- sample(non_tnbc_samples, length(valid_samples))

message("TNBC samples: ",       length(valid_samples))
message("Control samples: ",    length(control_samples))
message("Pan-BRCA samples: ",   length(all_unique_samples))

# ---------- 3. Gene extraction helper -----------------------------------------
get_gene_exp <- function(symbol, ensg, target_samples) {
  idx <- grep(ensg, rownames(tpm_matrix), value = TRUE)[1]
  if (is.na(idx))
    idx <- grep(paste0("^", symbol, "(\\.|$|\\|)"),
                rownames(tpm_matrix), value = TRUE)[1]
  if (is.na(idx)) {
    message("WARNING: Gene not found: ", symbol)
    return(rep(NA_real_, length(target_samples)))
  }
  as.numeric(unlist(tpm_matrix[idx, target_samples]))
}

# ---------- 4. Build data frames per group ------------------------------------
make_df <- function(target) {
  data.frame(
    P2RX4 = get_gene_exp("P2RX4", "ENSG00000258982", target),
    VIM   = get_gene_exp("VIM",   "ENSG00000026025", target),
    CDH1  = get_gene_exp("CDH1",  "ENSG00000039068", target)
  )
}

data_pan  <- make_df(all_unique_samples)
data_ctrl <- make_df(control_samples)
data_tnbc <- make_df(valid_samples)

# ---------- 5. Panels ---------------------------------------------------------
scatter_panel <- function(data, x, y, color, title, subtitle) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = color, alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", fill = color, alpha = 0.2) +
    stat_cor(method = "spearman",
             label.x.npc = "left", label.y.npc = "top", size = 4.5) +
    labs(title    = title,
         subtitle = subtitle,
         x = "P2RX4 Expression (TPM)",
         y = "VIM (Mesenchymal) Expression (TPM)") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", color = if (color == "#E67E22") "#D35400" else "black"))
}

p_pan  <- scatter_panel(data_pan,  "P2RX4", "VIM", "#95A5A6",
                        paste0("A. Pan-BRCA (n=", length(all_unique_samples), ")"),
                        "All subtypes mixed: No correlation")

p_ctrl <- scatter_panel(data_ctrl, "P2RX4", "VIM", "#3498DB",
                        paste0("B. Non-TNBC Control (n=", length(control_samples), ")"),
                        "Random sampled: No correlation")

p_tnbc <- scatter_panel(data_tnbc, "P2RX4", "VIM", "#E67E22",
                        paste0("C. TNBC Subtype (n=", length(valid_samples), ")"),
                        "Unique mechanism: Significant negative correlation!")

final_fig6 <- ggarrange(p_pan, p_ctrl, p_tnbc, ncol = 3, nrow = 1)

# ---------- 6. Save -----------------------------------------------------------
file_base <- file.path(OUTPUT_DIR, "Figure6_Subtype_Specificity")

ggsave(paste0(file_base, ".tiff"), plot = final_fig6,
       width = 15, height = 5, dpi = 300, compression = "lzw")
ggsave(paste0(file_base, ".pdf"),  plot = final_fig6,
       width = 15, height = 5)

message("✅ Figure 6 saved to: ", OUTPUT_DIR)

sessionInfo()

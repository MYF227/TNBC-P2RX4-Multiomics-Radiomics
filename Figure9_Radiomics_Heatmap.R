# =============================================================================
# Figure 9 — Heatmap of top peritumoral radiomics features
#             across 4 typical TNBC cases (3-mm peritumoral region)
#
# Data source : RADIOMICS_CSV (long-format, set in config.R)
#               Columns expected: CaseID, Feature_Class, Feature_Name,
#                                 GTV, Habitat_3mm
# Output      : Figure9_Radiomics_Heatmap.pdf
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(tidyverse)   # includes dplyr, tidyr (pivot_wider), readr
library(pheatmap)
library(RColorBrewer)

# ---------- 1. P2X4 TPM for the 4 pilot cases --------------------------------
p2x4_data <- data.frame(
  CaseID = c("A0E0", "A12F", "A0RX", "A0B3"),
  TPM    = c(2.68,   3.81,   8.31,   15.3)
) %>%
  mutate(Plot_Label = paste0(CaseID, "\n(", TPM, " TPM)"))

# ---------- 2. Load radiomics long-format data --------------------------------
if (!file.exists(RADIOMICS_CSV)) {
  stop("Radiomics CSV not found: ", RADIOMICS_CSV,
       "\nCheck RADIOMICS_CSV in config.R")
}

df <- read_csv(RADIOMICS_CSV, show_col_types = FALSE)

# ---------- 3. Compute Habitat / GTV ratio & clean labels --------------------
df_clean <- df %>%
  mutate(
    Ratio = Habitat_3mm / GTV,
    Clean_Class = case_when(
      Feature_Class == "firstorder" ~ "First-order",
      TRUE ~ toupper(Feature_Class)
    ),
    Clean_Feature_Name = paste0(Feature_Name, " (", Clean_Class, ")")
  ) %>%
  left_join(p2x4_data, by = "CaseID")

# ---------- 4. Wide format for correlation ------------------------------------
# pivot_wider replaces the deprecated spread()
wide_ratio <- df_clean %>%
  select(Plot_Label, Clean_Feature_Name, Ratio, TPM) %>%
  pivot_wider(names_from  = Clean_Feature_Name,
              values_from = Ratio)

# ---------- 5. Spearman correlation of each feature with P2X4 TPM -----------
feature_cols <- setdiff(names(wide_ratio), c("Plot_Label", "TPM"))

cor_results <- sapply(wide_ratio[, feature_cols], function(x) {
  cor(x, wide_ratio$TPM, method = "spearman", use = "complete.obs")
})

# Top 10 features (or fewer if less available)
n_top        <- min(10, length(cor_results))
top_features <- names(sort(abs(cor_results), decreasing = TRUE)[seq_len(n_top)])

message("Top ", n_top, " features selected for heatmap")

# ---------- 6. Build matrix sorted by P2X4 TPM (low → high) -----------------
wide_ratio <- wide_ratio %>% arrange(TPM)

plot_matrix            <- as.matrix(wide_ratio[, top_features])
rownames(plot_matrix)  <- wide_ratio$Plot_Label

# ---------- 7. Draw & save heatmap -------------------------------------------
out_pdf <- file.path(OUTPUT_DIR, "Figure9_Radiomics_Heatmap.pdf")

pheatmap(
  t(plot_matrix),
  scale        = "row",
  cluster_cols = FALSE,
  color        = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  main         = paste0(
    "P2X4 Expression vs Peritumoral Habitat Ratios\n",
    "(Pilot Analysis, n=4)"
  ),
  fontsize_row  = 10,
  fontsize_col  = 11,
  angle_col     = 0,
  border_color  = "white",
  filename      = out_pdf,
  width         = 8,
  height        = 6
)

message("✅ Figure 9 heatmap saved to: ", out_pdf)

sessionInfo()

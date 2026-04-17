# =============================================================================
# Figure 8A — Radiomics feature correlation with P2X4 expression
#              (Intratumoral GTV vs Peritumoral regions)
#
# Data source : Radiomics long-format CSV (RADIOMICS_CSV in config.R)
#               4 pilot cases: A0E0, A12F, A0RX, A0B3
#
# Prerequisite: GTV_df — a wide-format data frame of intratumoral GTV
#               radiomics features (cases × features), with a "CaseID" column.
#               Load or build this from your radiomics extraction pipeline
#               before running this script.
#
# Output      : Figure8A_Radiomics_Correlation scatter plot (PNG)
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(ggplot2)

# ---------- P2X4 TPM values for the 4 pilot cases ----------------------------
p2x4_values <- data.frame(
  CaseID   = c("A12F", "A0RX", "A0B3", "A0E0"),
  P2X4_TPM = c(3.81,   8.31,   15.3,   2.68)
)

# ---------- PREREQUISITE CHECK ------------------------------------------------
if (!exists("GTV_df")) {
  stop(
    "'GTV_df' not found in the environment.\n",
    "Please load your radiomics wide-format data frame first.\n",
    "It should have columns: CaseID, and one column per radiomics feature."
  )
}

# ---------- 1. Merge expression with radiomics --------------------------------
analysis_df <- merge(GTV_df, p2x4_values, by = "CaseID")

# ---------- 2. Identify numeric radiomics feature columns ---------------------
# Exclude CaseID and P2X4_TPM; keep only numeric columns
feature_cols <- names(analysis_df)[
  sapply(analysis_df, is.numeric) & names(analysis_df) != "P2X4_TPM"
]
message("Total numeric feature columns: ", length(feature_cols))

# ---------- 3. Spearman correlation with P2X4 ---------------------------------
cor_results <- sapply(analysis_df[, feature_cols, drop = FALSE], function(x) {
  cor(x, analysis_df$P2X4_TPM, method = "spearman", use = "complete.obs")
})

# Features with |rho| > 0.8
top_features <- names(cor_results[abs(cor_results) > 0.8])
message("Features with |rho| > 0.8: ", length(top_features))
if (length(top_features) == 0) {
  message("No features exceeded threshold; showing top 5 instead.")
  top_features <- names(sort(abs(cor_results), decreasing = TRUE))[1:5]
}

# ---------- 4. Example scatter: strongest correlated feature ------------------
best_feature <- top_features[1]

p_scatter <- ggplot(analysis_df,
                    aes(x = P2X4_TPM, y = .data[[best_feature]])) +
  geom_point(size = 4, color = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("P2X4 Expression vs.", best_feature),
    x     = "P2X4 mRNA Level (TPM)",
    y     = "Radiomics Feature Value"
  )

out_file <- file.path(OUTPUT_DIR, "Figure8A_Radiomics_Correlation.png")
ggsave(out_file, plot = p_scatter, width = 6, height = 5, dpi = 300)
message("✅ Figure 8A saved to: ", OUTPUT_DIR)

sessionInfo()

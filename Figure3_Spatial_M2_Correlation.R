# =============================================================================
# Figure 3 — Correlation between P2RX4 and M2 macrophage marker CD163
#
# Prerequisite : A Seurat spatial-transcriptomics object `st_obj` must be
#                loaded before running this script. Example:
#                  st_obj <- readRDS("data/spatial/CID44971_st_obj.rds")
#                or created from Visium raw data using Seurat's Load10X_Spatial().
#
# Output       : Figure3_P2RX4_CD163_Correlation.tiff / .eps
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(Seurat)
library(ggplot2)

# ---------- PREREQUISITE CHECK ------------------------------------------------
# ► Load your spatial transcriptomics Seurat object here before proceeding.
#   Replace the path below with the actual location of your saved object.

# st_obj <- readRDS("data/spatial/CID44971_st_obj.rds")
# -- OR, if building from raw Visium data: --
# st_obj <- Load10X_Spatial("data/spatial/CID44971_visium/")
# st_obj <- SCTransform(st_obj, assay = "Spatial", verbose = FALSE)

if (!exists("st_obj")) {
  stop(
    "Object 'st_obj' not found.\n",
    "Please load or build your spatial Seurat object before running this script.\n",
    "See the PREREQUISITE section at the top of this file for guidance."
  )
}

# ---------- 1. Extract expression data ----------------------------------------
# CD163  = classic M2 macrophage marker
# MRC1   = CD206, an additional M2 marker
exp_macro <- FetchData(st_obj, vars = c("P2RX4", "CD163", "MRC1"))

# ---------- 2. Spearman correlation -------------------------------------------
cor_m2 <- cor.test(exp_macro$P2RX4, exp_macro$CD163, method = "spearman")

cat("--- P2RX4 vs M2 Macrophage (CD163) Spatial Correlation ---\n")
print(cor_m2)

# ---------- 3. Scatter plot ---------------------------------------------------
p_cor <- ggplot(exp_macro, aes(x = P2RX4, y = CD163)) +
  geom_jitter(alpha = 0.3, color = "darkgreen") +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(
    title    = "TNBC CID44971: P2RX4 vs CD163 (M2 Macrophage)",
    subtitle = paste0(
      "Rho = ", round(cor_m2$estimate, 3),
      "  P = ", signif(cor_m2$p.value, 3)
    ),
    x = "P2RX4 Expression",
    y = "CD163 Expression"
  )

# ---------- 4. Save -----------------------------------------------------------
file_base <- file.path(OUTPUT_DIR, "Figure3_P2RX4_CD163_Correlation")

ggsave(paste0(file_base, ".tiff"), plot = p_cor,
       width = 7, height = 6, dpi = 300, compression = "lzw")
ggsave(paste0(file_base, ".eps"),  plot = p_cor,
       device = "eps", width = 7, height = 6)

message("✅ Figure 3 saved to: ", OUTPUT_DIR)

sessionInfo()

# =============================================================================
# Figure 5 — Significant enrichment of hypoxia-related pathways in
#             the P2RX4 high-expression group
#
# Prerequisites (must be created upstream, e.g. in a data-preparation script):
#   counts_clean  — normalised count matrix (genes × TNBC samples),
#                   extreme outliers already removed
#   group         — factor vector (levels: "Low", "High") corresponding to
#                   P2RX4 expression group for each sample in counts_clean
#   m_t2g_final   — MSigDB Hallmark gene-set data frame (TERM2GENE format),
#                   obtainable via:
#                     library(msigdbr)
#                     m_t2g_final <- msigdbr(species = "Homo sapiens",
#                                            category = "H")[, c("gs_name","gene_symbol")]
#
# Output : Figure5_GSEA_Hypoxia.png
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(limma)
library(clusterProfiler)   # provides GSEA()
library(enrichplot)        # provides gseaplot2()
library(msigdbr)           # for Hallmark gene sets (if not already loaded)

# ---------- PREREQUISITE CHECK ------------------------------------------------
required_objects <- c("counts_clean", "group", "m_t2g_final")
missing_objs     <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objs) > 0) {
  stop(
    "The following objects are required but not found in the environment:\n  ",
    paste(missing_objs, collapse = ", "), "\n",
    "Please run your upstream data-preparation script first.\n\n",
    "Quick guide to create m_t2g_final:\n",
    "  library(msigdbr)\n",
    "  m_t2g_final <- msigdbr(species = 'Homo sapiens', category = 'H')[,\n",
    "                   c('gs_name', 'gene_symbol')]"
  )
}

# ---------- 1. Differential expression (High vs Low) -------------------------
design      <- model.matrix(~ group)
fit         <- lmFit(log2(counts_clean + 1), design)
fit         <- eBayes(fit)
all_genes   <- topTable(fit, coef = 2, number = Inf)   # coef 2 = High - Low

# ---------- 2. Ranked gene list for GSEA ------------------------------------
gene_list <- sort(setNames(all_genes$logFC, rownames(all_genes)),
                  decreasing = TRUE)

# ---------- 3. Run GSEA -------------------------------------------------------
gsea_res <- GSEA(gene_list, TERM2GENE = m_t2g_final, pvalueCutoff = 0.05)

# ---------- 4. Plot HALLMARK_HYPOXIA enrichment plot -------------------------
pathway_id <- "HALLMARK_HYPOXIA"

if (pathway_id %in% gsea_res@result$ID) {

  p_gsea <- gseaplot2(
    gsea_res,
    geneSetID   = pathway_id,
    title       = "GSEA: Hallmark Hypoxia (High vs Low P2RX4)",
    pvalue_table = TRUE
  )

  out_file <- file.path(OUTPUT_DIR, "Figure5_GSEA_Hypoxia.png")
  ggsave(out_file, plot = p_gsea, width = 8, height = 6, dpi = 300)

  nes_val <- gsea_res@result[pathway_id, "NES"]
  message("✅ NES for HALLMARK_HYPOXIA: ", round(nes_val, 4))
  message("✅ Figure 5 saved to: ", OUTPUT_DIR)

} else {
  warning("HALLMARK_HYPOXIA not found in GSEA results at pvalueCutoff = 0.05.",
          "\nAvailable pathways:\n",
          paste(head(gsea_res@result$ID, 20), collapse = "\n"))
}

sessionInfo()

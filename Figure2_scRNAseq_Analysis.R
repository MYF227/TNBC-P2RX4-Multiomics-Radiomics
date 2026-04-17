# =============================================================================
# Figure 2 — Single-cell transcriptomic analysis of P2RX4 expression in TNBC
#
# Data source : GSE176078 (10x Genomics, CID44991 = TNBC, CID4535 = ER+)
# Output      : Figure2_scRNAseq_TNBC.tiff / .eps
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# ---------- Helper: load & process one 10x sample ----------------------------
fast_process <- function(mtx_path, sample_name) {
  message("Processing sample: ", sample_name)
  dir_now <- dirname(mtx_path)

  counts <- readMM(mtx_path)

  gene_f <- list.files(dir_now, pattern = "genes.tsv", full.names = TRUE)[1]
  cell_f <- list.files(dir_now, pattern = "barcodes.tsv", full.names = TRUE)[1]

  genes <- read.table(gene_f, stringsAsFactors = FALSE)
  cells <- read.table(cell_f, stringsAsFactors = FALSE)

  rownames(counts) <- make.unique(
    as.character(if (ncol(genes) >= 2) genes[, 2] else genes[, 1])
  )
  colnames(counts) <- as.character(cells[, 1])

  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj,          verbose = FALSE)
  obj <- FindVariableFeatures(obj,   selection.method = "vst",
                               nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj,              verbose = FALSE)
  obj <- RunPCA(obj,   npcs = 30,   verbose = FALSE)
  obj <- RunUMAP(obj,  dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  return(obj)
}

# ---------- 1. Process samples ------------------------------------------------
# ► Set correct sub-paths within SCRNA_DATA_DIR to match your folder layout
obj_tnbc <- fast_process(
  file.path(SCRNA_DATA_DIR, "GSM5354532_CID44991/CID44991/count_matrix_sparse.mtx"),
  "TNBC_44991"
)
obj_er <- fast_process(
  file.path(SCRNA_DATA_DIR, "GSM5354538_CID4535/CID4535/count_matrix_sparse.mtx"),
  "ER_4535"
)

# ---------- 2. Cell-type annotation -------------------------------------------
name_map <- c(
  "0"  = "Epithelial/Tumor", "1"  = "Epithelial/Tumor",
  "3"  = "Epithelial/Tumor", "6"  = "Epithelial/Tumor",
  "7"  = "Epithelial/Tumor", "13" = "Epithelial/Tumor",
  "15" = "Epithelial/Tumor",
  "4"  = "T cell",           "10" = "T cell",  "14" = "T cell",
  "11" = "Myeloid",
  "2"  = "B cell / Plasma cell", "5"  = "B cell / Plasma cell",
  "8"  = "B cell / Plasma cell", "12" = "B cell / Plasma cell",
  "9"  = "Fibroblast",
  "16" = "Endothelial"
)

rename_safe <- function(obj, mapping) {
  curr_lvls        <- levels(obj)
  matched          <- mapping[curr_lvls]
  matched[is.na(matched)] <- "Other"
  names(matched)   <- curr_lvls
  RenameIdents(obj, matched)
}

obj_tnbc_named <- rename_safe(obj_tnbc, name_map)
obj_er_named   <- rename_safe(obj_er,   name_map)

# ---------- 3. Panels ---------------------------------------------------------
gene <- "P2RX4"

p1 <- DimPlot(obj_tnbc_named, reduction = "umap", label = TRUE,
              label.size = 4, repel = TRUE) +
      theme_bw() + labs(title = "TNBC (Sample 44991)") + NoLegend()

p2 <- DimPlot(obj_er_named, reduction = "umap", label = TRUE,
              label.size = 4, repel = TRUE) +
      theme_bw() + labs(title = "ER+ (Sample 4535)") + NoLegend()

p3 <- FeaturePlot(obj_tnbc_named, features = gene,
                  cols = c("lightgrey", "red"), pt.size = 0.5) +
      theme_bw() + labs(title = paste(gene, "in TNBC"))

p4 <- FeaturePlot(obj_er_named, features = gene,
                  cols = c("lightgrey", "red"), pt.size = 0.5) +
      theme_bw() + labs(title = paste(gene, "in ER+"))

# ---------- 4. Compose & save ------------------------------------------------
fig2 <- (p1 | p2) / (p3 | p4) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 16, face = "bold"))

file_base <- file.path(OUTPUT_DIR, "Figure2_scRNAseq_TNBC")

ggsave(paste0(file_base, ".tiff"), plot = fig2,
       width = 14, height = 12, dpi = 300, compression = "lzw")
ggsave(paste0(file_base, ".eps"),  plot = fig2,
       device = "eps", width = 14, height = 12)

message("✅ Figure 2 saved to: ", OUTPUT_DIR)

# ---------- Session info (for reproducibility) --------------------------------
sessionInfo()

# =============================================================================
# Figure 7 / Figure 8B — Kaplan-Meier Overall Survival Analysis
#                         Pan-BRCA: P2RX4 High vs Low
#
# NOTE: The original document contained the complete survival analysis code
#       under "Figure 8B". It is presented here as the standalone Figure 7
#       script for clarity.
#
# Data source : TCGA-BRCA clinical TSV + TPM matrix (UCSC Xena / GDC portal)
# Output      : Figure7_P2RX4_PanBRCA_OS.tiff / .eps / .pdf (preview)
# =============================================================================

source("config.R")

# ---------- Packages ----------------------------------------------------------
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

# ---------- 1. Load clinical survival data ------------------------------------
cat("--- Step 1: Load clinical data ---\n")

cli_main  <- fread(CLINICAL_TSV, data.table = FALSE)

surv_data <- cli_main %>%
  mutate(
    CaseID = toupper(substr(submitter_id, 1, 12)),
    event  = ifelse(
      grepl("dead|deceased", vital_status.demographic, ignore.case = TRUE),
      1L, 0L
    ),
    time   = pmax(
      suppressWarnings(as.numeric(days_to_death.demographic)),
      suppressWarnings(as.numeric(days_to_last_follow_up.diagnoses)),
      na.rm = TRUE
    ) / 30.44      # convert days to months
  ) %>%
  filter(!is.na(time) & time > 0) %>%
  distinct(CaseID, .keep_all = TRUE)

cat("✅ Clinical records: ", nrow(surv_data), "\n")

# ---------- 2. Extract P2RX4 expression from TPM matrix ----------------------
cat("\n--- Step 2: Load TPM matrix (please wait) ---\n")

tpm_file <- paste0(TCGA_DATA_DIR, "/TCGA-BRCA.star_tpm.tsv")
if (!file.exists(tpm_file)) tpm_file <- sub("\\.tsv$", ".txt", tpm_file)
if (!file.exists(tpm_file)) stop("TPM file not found. Check TCGA_DATA_DIR in config.R")

exp_raw <- fread(tpm_file, data.table = FALSE)

# Locate P2RX4 row (search first 3 columns for gene symbol or Ensembl ID)
p2_idx <- which(apply(exp_raw[, 1:3], 1, function(x)
  any(grepl("P2RX4|ENSG00000100346", x, ignore.case = TRUE))))

if (length(p2_idx) == 0) stop("P2RX4 not found in TPM file. Check gene ID.")

p2_exp <- exp_raw[p2_idx[1], ] %>%
  pivot_longer(cols = -c(1:3), names_to = "SampleID", values_to = "Val") %>%
  mutate(
    SampleID_clean = gsub("\\.", "-", SampleID),
    SampleType     = substr(SampleID_clean, 14, 15),
    CaseID         = toupper(substr(SampleID_clean, 1, 12)),
    P2X4           = log2(as.numeric(Val) + 1)
  ) %>%
  filter(SampleType %in% paste0("0", 1:9)) %>%
  distinct(CaseID, .keep_all = TRUE)

cat("✅ P2RX4 expression extracted\n")

# ---------- 3. Merge ---------------------------------------------------------
cat("\n--- Step 3: Merge datasets ---\n")

final_df <- inner_join(surv_data, p2_exp, by = "CaseID")
cat("✅ Samples for analysis: ", nrow(final_df), "\n")

# ---------- 4. Optimal cutpoint & grouping ------------------------------------
res_cut   <- surv_cutpoint(final_df, time = "time", event = "event",
                            variables = "P2X4")
final_cat <- as.data.frame(surv_categorize(res_cut))
final_cat$P2X4 <- factor(final_cat$P2X4,
                          levels = c("high", "low"),
                          labels = c("P2X4_Value=high", "P2X4_Value=low"))

fit <- survfit(Surv(time, event) ~ P2X4, data = final_cat)
cat("Group sizes:\n"); print(table(final_cat$P2X4))

# ---------- 5. KM plot function -----------------------------------------------
draw_km <- function() {
  ggsurvplot(
    fit,
    data              = final_cat,
    pval              = TRUE,
    pval.size         = 5,
    pval.coord        = c(5, 0.18),
    conf.int          = TRUE,
    conf.int.alpha    = 0.15,
    palette           = c("#E7B800", "#2E9FDF"),
    risk.table        = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    tables.theme      = theme_cleantable(),
    title             = paste0("Pan-BRCA: P2RX4 OS (n=", nrow(final_cat), ")"),
    xlab              = "Time (Months)",
    ylab              = "Overall Survival Probability",
    legend.title      = "Strata",
    legend.labs       = c("P2X4_Value=high", "P2X4_Value=low"),
    font.main         = c(14, "bold"),
    font.x            = c(12, "plain"),
    font.y            = c(12, "plain"),
    font.tickslab     = c(10, "plain"),
    font.legend       = c(10, "plain"),
    break.time.by     = 50,
    xlim              = c(0, 300),
    censor.shape      = "|",
    censor.size       = 3,
    ggtheme = theme_classic(base_size = 12) +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.line       = element_line(linewidth = 0.6, color = "black"),
        axis.ticks      = element_line(linewidth = 0.5),
        panel.grid      = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.6, "cm"),
        plot.margin     = margin(10, 15, 5, 10)
      )
  )
}

# ---------- 6. Export TIFF (300 DPI) -----------------------------------------
out_tiff <- file.path(OUTPUT_DIR, "Figure7_P2RX4_PanBRCA_OS.tiff")

tiff(filename = out_tiff, width = 18, height = 20,
     units = "cm", res = 300, compression = "lzw")
print(draw_km())
dev.off()
cat("✅ TIFF saved: ", out_tiff, "\n")

# ---------- 7. Export EPS (vector) -------------------------------------------
out_eps <- file.path(OUTPUT_DIR, "Figure7_P2RX4_PanBRCA_OS.eps")

tryCatch({
  cairo_ps(filename = out_eps, width = 7.09, height = 7.87, onefile = FALSE)
  print(draw_km())
  dev.off()
  cat("✅ EPS saved: ", out_eps, "\n")
}, error = function(e) {
  postscript(file = out_eps, width = 7.09, height = 7.87,
             paper = "special", onefile = FALSE, horizontal = FALSE)
  print(draw_km())
  dev.off()
  cat("✅ EPS saved (postscript fallback): ", out_eps, "\n")
})

# ---------- 8. PDF preview ---------------------------------------------------
out_pdf <- file.path(OUTPUT_DIR, "Figure7_P2RX4_PanBRCA_OS_preview.pdf")
pdf(file = out_pdf, width = 7.09, height = 7.87)
print(draw_km())
dev.off()

cat("\n========== All done! ==========\n")
cat("Output folder: ", OUTPUT_DIR, "\n")

sessionInfo()

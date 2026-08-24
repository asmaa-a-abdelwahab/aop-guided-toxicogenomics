#!/usr/bin/env Rscript
# One RNA-seq dataset -> DESeq2 -> UP/DOWN/ALL gene sets and QC plots.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "matrixStats",
  "readr", "dplyr", "tidyr"
))
assert_files(c(COUNTS_FILE, META_FILE))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# --- load data ---------------------------------------------------------------
cts  <- read.delim(COUNTS_FILE, header = TRUE, row.names = 1, check.names = FALSE)
cts  <- cts[, vapply(cts, is.numeric, logical(1))]
meta <- readr::read_csv(META_FILE, show_col_types = FALSE)
stopifnot(all(c("sample_name","condition") %in% names(meta)))

rownames(meta) <- meta$sample_name
keep <- intersect(colnames(cts), meta$sample_name)
if (!length(keep)) {
  stop("No count-matrix columns match metadata sample_name values.", call. = FALSE)
}
cts  <- cts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
meta$condition <- dplyr::recode(meta$condition, "Medium"="CTRL", "Control"="CTRL", .default = meta$condition)
if (!"CTRL" %in% meta$condition) {
  stop("Metadata must contain a CTRL, Control, or Medium condition.", call. = FALSE)
}

# --- output dirs -------------------------------------------------------------
outdir <- DE_DIR
figdir <- file.path(DE_DIR, "figs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

# --- DESeq2 -------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(cts)), meta, ~ condition)
pref <- rowSums(DESeq2::counts(dds) >= 5) >= max(2, min(table(meta$condition)))
dds  <- dds[pref, ]
dds  <- DESeq2::DESeq(dds)

# Variance-stabilized matrix for downstream viz
vsd <- DESeq2::vst(dds, blind = TRUE)

# --- contrasts & exports ------------------------------------------------------
available_chems <- intersect(CHEMS, setdiff(unique(meta$condition), "CTRL"))
missing_chems <- setdiff(CHEMS, available_chems)
if (length(missing_chems)) {
  warning("Skipping conditions absent from metadata: ", paste(missing_chems, collapse = ", "))
}
if (!length(available_chems)) {
  stop("None of the configured CHEMS are present in metadata$condition.", call. = FALSE)
}
contrasts <- setNames(lapply(available_chems, \(x) c("condition", x, "CTRL")), available_chems)

# Collect per-contrast summary for barplot
deg_summaries <- list()

for (chem in names(contrasts)) {
  res <- DESeq2::results(dds, contrast = contrasts[[chem]])
  df  <- as.data.frame(res)
  df$ENSEMBL <- strip_v(rownames(df))

  # Save ALL
  readr::write_csv(df, file.path(outdir, sprintf("DEG_%s_vs_CTRL_all.csv", chem)))

  # Significant by thresholds
  sig <- subset(df, !is.na(padj) & padj < FDR_CUTOFF & abs(log2FoldChange) >= LFC_MIN)
  readr::write_csv(sig, file.path(outdir, sprintf("DEG_%s_vs_CTRL_significant.csv", chem)))

  # Lists for enrichment (choose UP/DOWN/ALL)
  if ("UP"   %in% RUN_ON)   readr::write_lines(unique(sig$ENSEMBL[sig$log2FoldChange > 0]),
                                               file.path(outdir, sprintf("DEG_%s_UP_ensg.txt", chem)))
  if ("DOWN" %in% RUN_ON)   readr::write_lines(unique(sig$ENSEMBL[sig$log2FoldChange < 0]),
                                               file.path(outdir, sprintf("DEG_%s_DOWN_ensg.txt", chem)))
  if ("ALL"  %in% RUN_ON)   readr::write_lines(unique(sig$ENSEMBL),
                                               file.path(outdir, sprintf("DEG_%s_ALL_ensg.txt", chem)))

  # -------- Volcano plot (per contrast) --------------------------------------
  # Use safe y transform and NA handling
  volcano_df <- df
  volcano_df$padj_plot <- ifelse(is.na(volcano_df$padj), 1, pmax(pmin(volcano_df$padj, 1), 1e-300))
  volcano_df$neglog10padj <- -log10(volcano_df$padj_plot)
  volcano_df$is_sig <- !is.na(volcano_df$padj) &
                       volcano_df$padj < FDR_CUTOFF &
                       abs(volcano_df$log2FoldChange) >= LFC_MIN

  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(aes(color = is_sig), size = 0.8, alpha = 0.8) +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70")) +
    theme_minimal(base_size = 11) +
    ggtitle(sprintf("Volcano: %s vs CTRL", chem)) +
    xlab("log2 Fold Change") + ylab("-log10 adjusted p-value") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(figdir, sprintf("Volcano_%s_vs_CTRL.png", chem)),
         p_volcano, width = 6, height = 5, dpi = 800)
  ggsave(file.path(figdir, sprintf("Volcano_%s_vs_CTRL.pdf", chem)),
         p_volcano, width = 6, height = 5)

  # -------- Add to summary ----------------------------------------------------
  n_up   <- sum(!is.na(sig$padj) & sig$log2FoldChange > 0)
  n_down <- sum(!is.na(sig$padj) & sig$log2FoldChange < 0)
  deg_summaries[[chem]] <- data.frame(
    Contrast = chem,
    n_up = n_up,
    n_down = n_down,
    n_sig = n_up + n_down,
    n_tested = sum(!is.na(df$padj)),
    stringsAsFactors = FALSE
  )
}

# ============================ Summary + Barplot ================================
deg_summary <- dplyr::bind_rows(deg_summaries)
alpha  <- FDR_CUTOFF
lfc_min <- LFC_MIN

# Save summary CSV
summary_file <- file.path(outdir, "DEG_summary_counts.csv")
write.csv(deg_summary, summary_file, row.names = FALSE)
print(deg_summary)

# Barplot (significant only)
deg_summary_long <- deg_summary |>
  tidyr::pivot_longer(cols = c("n_up", "n_down"),
                      names_to = "Direction", values_to = "Count")

p_bar <- ggplot(deg_summary_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("n_up" = "firebrick", "n_down" = "steelblue")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")) +
  ggtitle(sprintf("Significant DEGs per contrast (FDR < %.02f & |LFC| ≥ %s)", alpha, lfc_min))
ggsave(file.path(figdir, "DEG_summary_barplot.png"), p_bar, width = 7, height = 5, dpi = 800)
ggsave(file.path(figdir, "DEG_summary_barplot.pdf"), p_bar, width = 7, height = 5)

# ======================== Heatmap of Top Variable Genes ========================
# Use top 500 most variable genes from VST; center rows
topVarGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 500)
mat <- assay(vsd)[topVarGenes, , drop = FALSE]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
colnames(anno) <- "Condition"

# PNG
pheatmap(mat,
         annotation_col = anno,
         show_rownames = FALSE,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255),
         filename = file.path(figdir, "TopVarGenes_heatmap.png"),
         width = 8, height = 6)
# PDF
pdf(file.path(figdir, "TopVarGenes_heatmap.pdf"), width = 8, height = 6)
pheatmap(mat,
         annotation_col = anno,
         show_rownames = FALSE,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255))
dev.off()

cat("✓ DE + visualizations complete. Outputs in:\n  - Tables:", outdir, "\n  - Figures:", figdir, "\n")

# --- session info -------------------------------------------------------------
if (interactive()) {
  cat("\nSession info:\n")
  print(sessionInfo())
}

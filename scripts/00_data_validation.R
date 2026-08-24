#!/usr/bin/env Rscript
## ========================= RNA-seq Data Validation & QC =========================
## Outputs: outputs/Data_validation/
## ===============================================================================

## --- Setup ---------------------------------------------------------------------
source(file.path("scripts", "_config.R"))
assert_packages(c(
  "DESeq2", "data.table", "ggplot2", "pheatmap",
  "AnnotationDbi", "org.Hs.eg.db", "scales"
))
assert_files(c(COUNTS_FILE, META_FILE))

suppressPackageStartupMessages({
  library(data.table); library(DESeq2)
  library(ggplot2); library(pheatmap)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

## --- Paths & expected column names --------------------------------------------
OUT_DIR <- file.path(OUT_DIR, "Data_validation")

SAMPLE_COLS   <- c("sample_name","Run","run_accession","Sample","Sample.Name")
CONDITION_COL <- c("condition","group","Group","Condition")
BATCH_COLS    <- c("batch","Batch")  # optional

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## --- Provenance supplied by authors ------
upstream <- list(
  cellType = "PBMCs",
  time = "T24",
  rnaExtraction = "Tempus Spin RNA Isolation Kit (ThermoFisher)",
  contaminantsRemoval = "TurboDNAse kit (ThermoFisher)",
  rnaQuantification = "Qubit instrument (ThermoFisher)",
  rnaQC = "Agilent 2100 Bioanalyzer (Agilent, Santa Clara, California, USA)",
  mrnaSequencing = "QuantSeq 3’ mRNA-Seq Library Prep Kit FWD for Illumina (75 single-end)",
  
  instrument = "Illumina Hiseq4000",
  encoding      = "Sanger/Illumina 1.9",
  trimming      = "Adapters & polyA removed",
  aligner       = "STAR 2.6.1d",
  featureCounts = "v2.0.0",
  genomeAssembly = "GRCh38/hg38 (GCA_000001405.28, Dec 2013)"
)

## --- Helpers -------------------------------------------------------------------
is_integerish <- function(x) all(abs(x - round(x)) < 1e-8)
ggsave2 <- function(p, file, width=6, height=4, dpi=300) {
  suppressWarnings(ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi))
}

# DESeq2::plotPCA() currently triggers this ggplot2 lifecycle warning internally.
# Muffle only that known dependency warning; preserve all other warnings.
muffle_deseq2_aes_string_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      message <- conditionMessage(w)
      known_warning <-
        grepl("`aes_string()` was deprecated", message, fixed = TRUE) &&
        grepl("DESeq2", message, fixed = TRUE)
      if (known_warning) invokeRestart("muffleWarning")
    }
  )
}

write_lines <- function(x, file, append=TRUE) {
  cat(paste0(x, collapse=""), file=file, sep="", append=append)
}

QC_REPORT <- file.path(OUT_DIR, "QC_summary.md")
cat("# RNA-seq Data Validation Summary\n\n", file = QC_REPORT)

## --- Read counts ---------------------------------------------------------------
cts_dt <- fread(COUNTS_FILE, sep="\t", header=TRUE, check.names=FALSE)
stopifnot(ncol(cts_dt) >= 2)

gene_col <- names(cts_dt)[1]
genes    <- cts_dt[[gene_col]]

# Keep only numeric sample columns
num_cols <- names(cts_dt)[sapply(cts_dt, is.numeric)]
if (length(num_cols) == 0) stop("No numeric columns detected in counts file.")
cts <- as.matrix(cts_dt[, ..num_cols])
rownames(cts) <- genes

# Collapse duplicate gene IDs
n_dups <- sum(duplicated(rownames(cts)))
if (n_dups > 0) cts <- rowsum(cts, group=rownames(cts))

write_lines(sprintf("**Counts matrix:** %d genes × %d samples\n\n", nrow(cts), ncol(cts)), QC_REPORT)
write_lines(sprintf("- Duplicated gene IDs collapsed: %d\n\n", n_dups), QC_REPORT)

# Raw integer check
all_integer <- all(apply(cts, 2, is_integerish))
write_lines(sprintf("- All columns integer-ish (raw counts): %s\n\n", ifelse(all_integer,"YES","NO (warning)")), QC_REPORT)
if (!all_integer) warning("Counts are not all integers. These may have been normalized already.")

# Library sizes (raw)
lib_sizes_raw <- colSums(cts)
write.csv(data.frame(sample=names(lib_sizes_raw), lib_size=lib_sizes_raw),
          file.path(OUT_DIR, "library_sizes_raw.csv"), row.names = FALSE)
write_lines("**Library sizes (raw counts)** saved: `library_sizes_raw.csv`\n\n", QC_REPORT)
write_lines(c("```\n", capture.output(summary(lib_sizes_raw)), "\n```\n\n"), QC_REPORT)

## --- Read metadata -------------------------------------------------------------
meta0 <- fread(META_FILE)

sample_col <- SAMPLE_COLS[SAMPLE_COLS %in% names(meta0)][1]
if (is.na(sample_col)) stop("Could not find a sample ID column in metadata (tried: ",
                            paste(SAMPLE_COLS, collapse=", "), ").")
cond_col   <- CONDITION_COL[CONDITION_COL %in% names(meta0)][1]
if (is.na(cond_col)) stop("Could not find a condition column in metadata (tried: ",
                          paste(CONDITION_COL, collapse=", "), ").")
batch_col  <- BATCH_COLS[BATCH_COLS %in% names(meta0)][1]  # optional

meta <- as.data.frame(meta0[, c(sample_col, cond_col, batch_col)[c(TRUE, TRUE, !is.na(batch_col))], with=FALSE])
names(meta)[1:2] <- c("sample","condition")
if (!is.na(batch_col)) names(meta)[3] <- "batch"

# Ensure counts columns all exist in metadata
if (!all(colnames(cts) %in% meta$sample)) {
  missing <- setdiff(colnames(cts), meta$sample)
  stop("These count columns aren’t in metadata$sample: ", paste(missing, collapse=", "))
}
meta <- meta[match(colnames(cts), meta$sample), , drop=FALSE]
rownames(meta) <- meta$sample

# Drop samples with missing condition
keep_samp <- !is.na(meta$condition) & meta$condition != ""
cts  <- cts[, keep_samp, drop=FALSE]
meta <- meta[keep_samp, , drop=FALSE]

# Factorize design columns
meta$condition <- factor(meta$condition)
has_batch <- "batch" %in% names(meta) && any(!is.na(meta$batch))
if (has_batch) meta$batch <- factor(meta$batch)
design_formula <- as.formula(if (has_batch) "~ batch + condition" else "~ condition")

write_lines(sprintf("**Metadata rows after matching:** %d\n\n", nrow(meta)), QC_REPORT)
write.csv(meta, file.path(OUT_DIR, "metadata_matched.csv"), row.names=FALSE)

## --- Gene ID handling & mapping ------------------------------------------------
# Strip Ensembl version suffix and detect Ensembl-like IDs
rownames(cts) <- sub("\\.\\d+$","", rownames(cts))
is_ens <- grepl("^ENSG\\d+$", rownames(cts))
write_lines(sprintf("**Rows looking like Ensembl IDs:** %d / %d\n\n", sum(is_ens), nrow(cts)), QC_REPORT)

# Map to symbols (best-effort; not all Ensembl have SYMBOLs)
gene_symbol <- setNames(rep(NA_character_, nrow(cts)), rownames(cts))
if (any(is_ens)) {
  ens_ids <- unique(rownames(cts)[is_ens])
  map <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys=ens_ids, keytype="ENSEMBL", columns="SYMBOL")
  )
  map <- map[!duplicated(map$ENSEMBL), c("ENSEMBL","SYMBOL")]
  gene_symbol[map$ENSEMBL] <- map$SYMBOL
} else {
  gene_symbol <- rownames(cts)
}
sym <- unname(gene_symbol)
write_lines(sprintf("**Mapped gene symbols available for:** %d genes\n\n", sum(!is.na(sym))), QC_REPORT)
write.csv(data.frame(ENSEMBL=names(gene_symbol), SYMBOL=gene_symbol),
          file.path(OUT_DIR,"ensembl_to_symbol_map.csv"), row.names=FALSE)

## --- Filtering (CPM rule) -----------------------------------------------------
group <- factor(meta$condition)
min_reps <- max(2, min(table(group)))  # ≥2 by default
lib_sizes <- colSums(cts)
cpm_matrix <- sweep(cts, 2, lib_sizes / 1e6, "/")
keep <- rowSums(cpm_matrix > 1) >= min_reps
write_lines(sprintf("**CPM filter:** kept %d / %d genes (CPM>1 in ≥ %d samples)\n\n",
                    sum(keep), nrow(cts), min_reps), QC_REPORT)

## --- DESeq2 normalization / VST for QC ----------------------------------------
cts_filt <- cts[keep, , drop = FALSE]
storage.mode(cts_filt) <- "integer"
stopifnot(all(colnames(cts_filt) == rownames(meta)))
cts_filt <- cts_filt[rowSums(cts_filt) > 0, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = cts_filt, colData = meta, design = design_formula)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

write.csv(norm_counts, file.path(OUT_DIR, "normalized_counts_deseq2.csv"))
saveRDS(vsd_mat, file = file.path(OUT_DIR, "vst_matrix.rds"))
write_lines("**DESeq2:** size-factor normalization + dispersion fit completed.\n\n", QC_REPORT)

## --- Exploratory QC: PCA & correlations ---------------------------------------
# PCA
p_pca <- suppressMessages(muffle_deseq2_aes_string_warning(
  plotPCA(vsd, intgroup = intersect(c("condition","batch"), colnames(colData(vsd))))
)) +
  ggtitle("PCA (VST, top 500 HVGs)")
ggsave2(p_pca, file.path(OUT_DIR, "pca_vst.png"), width=6.5, height=5)

# Sample correlations
cors <- cor(vsd_mat)
png(file.path(OUT_DIR, "sample_correlations_vst.png"), width=950, height=850)
pheatmap(cors, main="Sample correlations (VST)", fontsize=10)
invisible(dev.off())

write_lines("- PCA plot saved: `pca_vst.png`\n", QC_REPORT)
write_lines("- Sample correlation heatmap saved: `sample_correlations_vst.png`\n\n", QC_REPORT)

## --- Biological sanity checks --------------------------------------------------
# 1) Mito / ribo fractions from **raw** counts (unfiltered)
cm0 <- cts                           # raw counts after metadata matching
ids0 <- sub("\\.\\d+$", "", rownames(cm0))
map0 <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = ids0, keytype = "ENSEMBL", columns = "SYMBOL")
)
map0 <- map0[!duplicated(map0$ENSEMBL), ]
sym0 <- map0$SYMBOL[match(ids0, map0$ENSEMBL)]

mt_idx0 <- !is.na(sym0) & grepl("^MT-", sym0)
rp_idx0 <- !is.na(sym0) & grepl("^(RPL|RPS)", sym0)
safe_colsum <- function(m, idx) if (!any(idx)) rep(0, ncol(m)) else colSums(m[idx,,drop=FALSE])

fractions0 <- rbind(
  mito = safe_colsum(cm0, mt_idx0) / colSums(cm0),
  ribo = safe_colsum(cm0, rp_idx0) / colSums(cm0)
)
write.csv(round(t(fractions0),4), file.path(OUT_DIR,"mito_ribo_fractions_from_raw.csv"))
write_lines("**Mito/Ribo fractions (from raw counts)** saved: `mito_ribo_fractions_from_raw.csv`\n\n", QC_REPORT)

# Quick barplots
fr_long <- data.frame(
  sample = rep(colnames(cm0), each=2),
  fraction = as.numeric(fractions0),
  type = rep(c("mito","ribo"), times = ncol(cm0))
)
p_frac <- ggplot(fr_long, aes(x = sample, y = fraction, fill = type)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title = "Mitochondrial / Ribosomal fractions (raw counts)", x="", y="fraction")
ggsave2(p_frac, file.path(OUT_DIR,"mito_ribo_barplot.png"), width=7.5, height=6)

# 2) Sex markers (normalized counts for readability)
wanted <- c("XIST","DDX3Y","UTY","RPS4Y1","KDM5D")
sym2rows <- split(seq_len(nrow(norm_counts)),
                  map0$SYMBOL[match(rownames(norm_counts), map0$ENSEMBL)])
getCounts <- function(g) {
  idx <- sym2rows[[g]]
  if (is.null(idx)) rep(NA_real_, ncol(norm_counts)) else colSums(norm_counts[idx, , drop=FALSE])
}
sex_counts <- sapply(wanted, getCounts, simplify="matrix")
rownames(sex_counts) <- colnames(norm_counts)
write.csv(sex_counts, file.path(OUT_DIR, "sex_marker_normalizedCounts.csv"))
write_lines("**Sex markers (normalized counts)** saved: `sex_marker_normalizedCounts.csv`\n\n", QC_REPORT)

## --- Light-touch outlier flags -------------------------------------------------
# Simple flags (tweak thresholds to taste)
flags <- data.frame(
  sample = colnames(cts),
  lib_size_raw = lib_sizes_raw[colnames(cts)],
  mito_frac = as.numeric(fractions0["mito", colnames(cts)]),
  ribo_frac = as.numeric(fractions0["ribo", colnames(cts)])
)
flags$low_depth  <- flags$lib_size_raw < quantile(flags$lib_size_raw, 0.10)  # bottom 10%
flags$high_mito  <- flags$mito_frac > 0.20                                   # >20%
flags$high_ribo  <- flags$ribo_frac > 0.20                                   # usually rare in bulk
write.csv(flags, file.path(OUT_DIR,"sample_flags.csv"), row.names=FALSE)
write_lines("**Sample flags** (low depth / high mito / high ribo) saved: `sample_flags.csv`\n\n", QC_REPORT)

## --- Provenance & session info -------------------------------------------------
# Write provenance & design
write_lines("## Provenance (as provided by authors)\n\n", QC_REPORT)
write_lines(paste0("- Cell type: ", upstream$cellType, "\n"), QC_REPORT)
write_lines(paste0("- Time point: ", upstream$time, "\n"), QC_REPORT)
write_lines(paste0("- RNA extraction: ", upstream$rnaExtraction, "\n"), QC_REPORT)
write_lines(paste0("- Contaminants removal: ", upstream$contaminantsRemoval, "\n"), QC_REPORT)
write_lines(paste0("- RNA quantification: ", upstream$rnaQuantification, "\n"), QC_REPORT)
write_lines(paste0("- RNA QC: ", upstream$rnaQC, "\n"), QC_REPORT)
write_lines(paste0("- mRNA sequencing: ", upstream$mrnaSequencing, "\n"), QC_REPORT)
write_lines(paste0("- Instrument: ", upstream$instrument, "\n"), QC_REPORT)
write_lines(paste0("- Encoding: ", upstream$encoding, "\n"), QC_REPORT)
write_lines(paste0("- Trimming: ", upstream$trimming, "\n"), QC_REPORT)
write_lines(paste0("- STAR: ", upstream$aligner, "\n"), QC_REPORT)
write_lines(paste0("- featureCounts: ", upstream$featureCounts, "\n"), QC_REPORT)
write_lines(paste0("- Genome assembly: ", upstream$genomeAssembly, "\n"), QC_REPORT)

write_lines("## Analysis design\n\n", QC_REPORT)
write_lines(paste0("- Design formula: `", deparse(design_formula), "`\n"), QC_REPORT)
write_lines(paste0("- Conditions & replicate counts:\n\n```\n",
                   paste(capture.output(table(meta$condition)), collapse="\n"),
                   "\n```\n\n"), QC_REPORT)
if (has_batch) {
  write_lines(paste0("- Batches:\n\n```\n",
                     paste(capture.output(table(meta$batch)), collapse="\n"),
                     "\n```\n\n"), QC_REPORT)
}

# Session info
sink(file.path(OUT_DIR, "sessionInfo.txt")); print(sessionInfo()); sink()
write_lines("**Session info** saved: `sessionInfo.txt`\n\n", QC_REPORT)

write_lines("---\n\n✅ Data validation completed. See `outputs/Data_validation/`.\n", QC_REPORT)


# ---- Library-size bar plot (like the attached figure) -------------------------
library(scales)

# Data (ordered as in counts/metadata)
samples <- colnames(cts)
df_lib  <- data.frame(
  sample        = samples,
  lib_millions  = as.numeric(lib_sizes_raw[samples]) / 1e6,
  stringsAsFactors = FALSE
)

# Optional compact x labels "S1..S18" and a mapping table (saved for reference)
df_lib$label <- paste0("S", seq_len(nrow(df_lib)))
write.csv(df_lib[, c("label","sample","lib_millions")],
          file.path(OUT_DIR, "library_sizes_millions_with_labels.csv"),
          row.names = FALSE)

# Depth threshold (in millions); change if you want a different line
depth_thr_m <- 2.0

p_lib <- ggplot(df_lib, aes(x = factor(label, levels = label), y = lib_millions)) +
  geom_col(fill = "#90CAF9", width = 0.82) +
  geom_hline(aes(yintercept = depth_thr_m, color = "Depth threshold"),
             linetype = "dashed", linewidth = 0.8, show.legend = TRUE) +
  scale_color_manual("", values = c("Depth threshold" = "red")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Library Sizes per Sample",
    x = NULL,
    y = "Library size (millions)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.16, 0.90),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Save (PNG + PDF)
ggsave2(p_lib, file.path(OUT_DIR, "library_sizes_bar.png"),  width = 10, height = 5, dpi = 300)
ggsave2(p_lib, file.path(OUT_DIR, "library_sizes_bar.pdf"),  width = 10, height = 5, dpi = 300)

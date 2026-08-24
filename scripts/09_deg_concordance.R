#!/usr/bin/env Rscript
# Cross-contrast DEG concordance, recurrence, and effect-size synthesis.
# Discovers step-01 outputs dynamically so it can be reused across case studies.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "AnnotationDbi", "org.Hs.eg.db", "readr", "dplyr", "tidyr", "tibble",
  "ggplot2", "pheatmap", "scales", "stringr"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

RESULT_DIR <- file.path(OUT_DIR, "Interpretation", "DEG_concordance")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

PUBLICATION_DPI <- suppressWarnings(as.integer(
  Sys.getenv("AOP_PUBLICATION_DPI", unset = "600")
))
if (!is.finite(PUBLICATION_DPI) || PUBLICATION_DPI < 72L) PUBLICATION_DPI <- 600L
TOP_GENE_COUNT <- suppressWarnings(as.integer(
  Sys.getenv("AOP_PUBLICATION_TOP_GENES", unset = "50")
))
if (!is.finite(TOP_GENE_COUNT) || TOP_GENE_COUNT < 1L) TOP_GENE_COUNT <- 50L

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(RESULT_DIR, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = PUBLICATION_DPI,
    bg = "white"
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    file.path(RESULT_DIR, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = pdf_device,
    bg = "white"
  )
}

all_files <- list.files(
  DE_DIR,
  pattern = "^DEG_.*_vs_CTRL_all\\.csv$",
  full.names = TRUE
)
if (!length(all_files)) {
  stop("No step-01 DE result tables found under: ", DE_DIR, call. = FALSE)
}

experiment_from_file <- function(path) {
  sub("^DEG_(.*)_vs_CTRL_all\\.csv$", "\\1", basename(path))
}

read_de <- function(path) {
  data <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  required <- c("ENSEMBL", "log2FoldChange", "pvalue", "padj")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(basename(path), " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data %>%
    transmute(
      Experiment = experiment_from_file(path),
      ENSEMBL = strip_v(as.character(ENSEMBL)),
      log2FoldChange = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj),
      significant = !is.na(padj) & padj < FDR_CUTOFF &
        is.finite(log2FoldChange) & abs(log2FoldChange) >= LFC_MIN,
      Direction = case_when(
        significant & log2FoldChange > 0 ~ "UP",
        significant & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "Not significant"
      )
    ) %>%
    filter(!is.na(ENSEMBL), nzchar(ENSEMBL))
}

de_long <- bind_rows(lapply(all_files, read_de))
experiment_order <- sort(unique(de_long$Experiment))
if (length(experiment_order) < 2L) {
  stop("Step 09 requires at least two DE contrasts for concordance analysis.", call. = FALSE)
}

valid_ensembl <- intersect(
  unique(de_long$ENSEMBL),
  AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
)
symbol_table <- if (length(valid_ensembl)) {
  suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_ensembl,
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  )) %>%
    filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
    arrange(ENSEMBL, SYMBOL) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
} else {
  data.frame(ENSEMBL = character(), SYMBOL = character())
}
symbol_map <- setNames(symbol_table$SYMBOL, symbol_table$ENSEMBL)
de_long$SYMBOL <- unname(symbol_map[de_long$ENSEMBL])

# ---- DEG burden -------------------------------------------------------------
deg_counts <- de_long %>%
  group_by(Experiment) %>%
  summarise(
    n_up = sum(Direction == "UP"),
    n_down = sum(Direction == "DOWN"),
    n_significant = n_up + n_down,
    n_tested = sum(!is.na(padj)),
    significant_fraction = n_significant / pmax(n_tested, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_significant))
readr::write_csv(deg_counts, file.path(RESULT_DIR, "deg_burden_by_experiment.csv"))

burden_plot_data <- deg_counts %>%
  select(Experiment, n_up, n_down) %>%
  pivot_longer(c(n_up, n_down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Direction = recode(Direction, n_up = "UP", n_down = "DOWN"),
    Signed_count = ifelse(Direction == "DOWN", -Count, Count),
    Experiment = factor(Experiment, levels = rev(deg_counts$Experiment))
  )

p_burden <- ggplot(burden_plot_data, aes(Experiment, Signed_count, fill = Direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = c(UP = "#B2182B", DOWN = "#2166AC")) +
  scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
  labs(
    title = "Differential-expression burden by contrast",
    subtitle = sprintf("FDR < %.2g and |log2 fold change| >= %s", FDR_CUTOFF, LFC_MIN),
    x = NULL,
    y = "Significant genes",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )
save_plot(p_burden, "Figure_DEG_burden", width = 8, height = 5.5)

# ---- Pairwise significant-gene overlap -------------------------------------
make_gene_sets <- function(direction) {
  filtered <- de_long %>% filter(significant)
  if (direction != "ALL") filtered <- filtered %>% filter(Direction == direction)
  split(filtered$ENSEMBL, filtered$Experiment) %>% lapply(unique)
}

pairwise_overlap <- function(gene_sets, direction) {
  experiments <- experiment_order
  tidyr::expand_grid(Experiment_A = experiments, Experiment_B = experiments) %>%
    rowwise() %>%
    mutate(
      n_A = length(gene_sets[[Experiment_A]] %||% character()),
      n_B = length(gene_sets[[Experiment_B]] %||% character()),
      n_intersection = length(intersect(
        gene_sets[[Experiment_A]] %||% character(),
        gene_sets[[Experiment_B]] %||% character()
      )),
      n_union = length(union(
        gene_sets[[Experiment_A]] %||% character(),
        gene_sets[[Experiment_B]] %||% character()
      )),
      Jaccard = ifelse(n_union > 0, n_intersection / n_union, NA_real_),
      Direction = direction
    ) %>%
    ungroup()
}

`%||%` <- function(x, y) if (is.null(x)) y else x
overlap_long <- bind_rows(lapply(c("ALL", "UP", "DOWN"), function(direction) {
  pairwise_overlap(make_gene_sets(direction), direction)
}))
readr::write_csv(overlap_long, file.path(RESULT_DIR, "pairwise_DEG_overlap.csv"))

overlap_plot_data <- overlap_long %>%
  mutate(
    Experiment_A = factor(Experiment_A, levels = experiment_order),
    Experiment_B = factor(Experiment_B, levels = rev(experiment_order)),
    Direction = factor(Direction, levels = c("ALL", "UP", "DOWN")),
    Label = ifelse(is.na(Jaccard), "NA", sprintf("%.2f", Jaccard))
  )

p_overlap <- ggplot(overlap_plot_data, aes(Experiment_A, Experiment_B, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 3) +
  facet_wrap(~Direction, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#762A83", limits = c(0, 1), na.value = "grey90") +
  labs(
    title = "Pairwise overlap of significant gene sets",
    subtitle = "Jaccard index; direction-specific panels separate concordant responses",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
save_plot(p_overlap, "Figure_DEG_Jaccard_overlap", width = 14, height = 5.2)

# ---- Genome-wide effect-size concordance -----------------------------------
lfc_wide <- de_long %>%
  select(ENSEMBL, Experiment, log2FoldChange) %>%
  distinct() %>%
  pivot_wider(names_from = Experiment, values_from = log2FoldChange)
lfc_matrix <- as.matrix(lfc_wide[, setdiff(names(lfc_wide), "ENSEMBL"), drop = FALSE])
correlation_matrix <- stats::cor(
  lfc_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)
readr::write_csv(
  as.data.frame(correlation_matrix) %>% tibble::rownames_to_column("Experiment"),
  file.path(RESULT_DIR, "log2FC_spearman_correlation.csv")
)

correlation_colors <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
pheatmap::pheatmap(
  correlation_matrix,
  color = correlation_colors,
  breaks = seq(-1, 1, length.out = 102),
  border_color = "white",
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Genome-wide log2FC concordance (Spearman rho)",
  filename = file.path(RESULT_DIR, "Figure_log2FC_correlation.png"),
  width = 7,
  height = 6
)
pheatmap::pheatmap(
  correlation_matrix,
  color = correlation_colors,
  breaks = seq(-1, 1, length.out = 102),
  border_color = "white",
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Genome-wide log2FC concordance (Spearman rho)",
  filename = file.path(RESULT_DIR, "Figure_log2FC_correlation.pdf"),
  width = 7,
  height = 6
)

# ---- Recurrent significant genes -------------------------------------------
recurrent_genes <- de_long %>%
  filter(significant) %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL[!is.na(SYMBOL)], default = NA_character_),
    n_experiments = n_distinct(Experiment),
    n_up = sum(Direction == "UP"),
    n_down = sum(Direction == "DOWN"),
    direction_consistency = pmax(n_up, n_down) / n_experiments,
    median_abs_log2FC = median(abs(log2FoldChange), na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    experiments = paste(sort(unique(Experiment)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_experiments), desc(direction_consistency), desc(median_abs_log2FC))
if (!nrow(recurrent_genes)) {
  stop(
    "No genes pass the configured FDR and fold-change thresholds; recurrence cannot be calculated.",
    call. = FALSE
  )
}
readr::write_csv(recurrent_genes, file.path(RESULT_DIR, "recurrent_significant_genes.csv"))

recurrence_distribution <- recurrent_genes %>%
  count(n_experiments, name = "n_genes") %>%
  arrange(n_experiments)
readr::write_csv(
  recurrence_distribution,
  file.path(RESULT_DIR, "gene_recurrence_distribution.csv")
)

p_recurrence <- ggplot(recurrence_distribution, aes(n_experiments, n_genes)) +
  geom_col(fill = "#4D4D4D", width = 0.72) +
  geom_text(aes(label = scales::comma(n_genes)), vjust = -0.25, size = 3.4) +
  scale_x_continuous(breaks = recurrence_distribution$n_experiments) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Recurrence of significant genes across contrasts",
    x = "Number of contrasts",
    y = "Genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold"))
save_plot(p_recurrence, "Figure_gene_recurrence", width = 6.5, height = 4.8)

top_gene_count <- min(TOP_GENE_COUNT, nrow(recurrent_genes))
if (top_gene_count > 1L) {
  top_genes <- recurrent_genes %>% slice_head(n = top_gene_count)
  gene_labels <- ifelse(
    is.na(top_genes$SYMBOL),
    top_genes$ENSEMBL,
    paste0(top_genes$SYMBOL, " | ", top_genes$ENSEMBL)
  )
  names(gene_labels) <- top_genes$ENSEMBL

  top_lfc <- de_long %>%
    filter(ENSEMBL %in% top_genes$ENSEMBL) %>%
    select(ENSEMBL, Experiment, log2FoldChange) %>%
    distinct() %>%
    pivot_wider(names_from = Experiment, values_from = log2FoldChange) %>%
    tibble::column_to_rownames("ENSEMBL") %>%
    as.matrix()
  rownames(top_lfc) <- unname(gene_labels[rownames(top_lfc)])
  top_lfc <- top_lfc[unname(gene_labels[top_genes$ENSEMBL]), , drop = FALSE]

  limit <- max(abs(top_lfc), na.rm = TRUE)
  pheatmap::pheatmap(
    top_lfc,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = seq(-limit, limit, length.out = 102),
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    main = sprintf("Top %d recurrent DEGs: log2 fold change", top_gene_count),
    filename = file.path(RESULT_DIR, "Figure_recurrent_gene_log2FC.png"),
    width = 9,
    height = 12
  )
  pheatmap::pheatmap(
    top_lfc,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = seq(-limit, limit, length.out = 102),
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    main = sprintf("Top %d recurrent DEGs: log2 fold change", top_gene_count),
    filename = file.path(RESULT_DIR, "Figure_recurrent_gene_log2FC.pdf"),
    width = 9,
    height = 12
  )
}

# ---- Data-driven interpretation summary ------------------------------------
off_diagonal <- correlation_matrix
diag(off_diagonal) <- NA_real_
best_pair <- which(off_diagonal == max(off_diagonal, na.rm = TRUE), arr.ind = TRUE)[1, ]
most_shared <- recurrent_genes %>% filter(n_experiments == max(n_experiments))

summary_lines <- c(
  "# DEG concordance summary",
  "",
  sprintf("- Contrasts analyzed: %d (%s).", length(experiment_order), paste(experiment_order, collapse = ", ")),
  sprintf(
    "- Largest DEG burden: %s (%s significant genes).",
    deg_counts$Experiment[[1]],
    scales::comma(deg_counts$n_significant[[1]])
  ),
  sprintf(
    "- Strongest genome-wide effect concordance: %s vs %s (Spearman rho = %.3f).",
    rownames(off_diagonal)[best_pair[[1]]],
    colnames(off_diagonal)[best_pair[[2]]],
    off_diagonal[best_pair[[1]], best_pair[[2]]]
  ),
  sprintf(
    "- Maximum DEG recurrence: %d contrasts (%d genes reach this recurrence).",
    max(recurrent_genes$n_experiments),
    nrow(most_shared)
  ),
  "",
  "These summaries prioritize recurrent and concordant signals; they do not establish causality."
)
writeLines(summary_lines, file.path(RESULT_DIR, "interpretation_summary.md"), useBytes = TRUE)

message("DEG concordance analysis complete: ", normalizePath(RESULT_DIR, winslash = "/"))

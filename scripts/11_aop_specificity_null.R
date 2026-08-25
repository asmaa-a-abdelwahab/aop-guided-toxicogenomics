#!/usr/bin/env Rscript
# Size-matched random-background analysis for KE/AOP enrichment specificity.

source(file.path("scripts", "_config.R"))
assert_packages(c("AOPfingerprintR", "Matrix", "readr", "readxl", "dplyr", "ggplot2", "scales"))

suppressPackageStartupMessages({
  library(AOPfingerprintR)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

NULL_ITERATIONS <- as.integer(numeric_setting("AOP_NULL_ITERATIONS", 1000))
NULL_SEED <- as.integer(numeric_setting("AOP_NULL_SEED", 20260825))
NULL_DIRECTIONS <- toupper(csv_setting("AOP_NULL_DIRECTIONS", "ALL"))
NULL_BATCH_SIZE <- as.integer(numeric_setting("AOP_NULL_BATCH_SIZE", 100))

if (NULL_ITERATIONS < 99L) stop("AOP_NULL_ITERATIONS must be at least 99.", call. = FALSE)
if (NULL_BATCH_SIZE < 1L) stop("AOP_NULL_BATCH_SIZE must be positive.", call. = FALSE)
if (!all(NULL_DIRECTIONS %in% c("ALL", "UP", "DOWN"))) {
  stop("AOP_NULL_DIRECTIONS may contain only ALL, UP, and DOWN.", call. = FALSE)
}

SPEC_DIR <- file.path(OUT_DIR, "Specificity_null")
dir.create(SPEC_DIR, recursive = TRUE, showWarnings = FALSE)

normalize_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  unique(x[!is.na(x) & nzchar(x)])
}

read_background <- function(path) {
  normalize_ids(suppressMessages(readxl::read_excel(path, col_names = FALSE))[[1]])
}

assert_files(c(BACKGROUND_GSE_FILE, BACKGROUND_LIT_FILE), "Background gene files")
backgrounds <- list(
  Transcriptomic = read_background(BACKGROUND_GSE_FILE),
  Literature = read_background(BACKGROUND_LIT_FILE)
)

read_gene_column <- function(path) {
  dat <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  candidates <- intersect(
    c("Gene", "Feature", "ENSEMBL", "Ensembl", "gene", "GeneID", "ID"),
    names(dat)
  )
  if (!length(candidates)) {
    stop("No gene column found in ", path, ".", call. = FALSE)
  }
  normalize_ids(dat[[candidates[[1]]]])
}

discover_gene_sets <- function() {
  records <- list()

  for (chem in CHEMS) {
    for (direction in intersect(RUN_ON, NULL_DIRECTIONS)) {
      path <- file.path(DE_DIR, sprintf("DEG_%s_%s_ensg.txt", chem, direction))
      if (!file.exists(path)) next
      key <- paste("Transcriptomic", chem, direction, sep = "::")
      records[[key]] <- list(
        branch = "Transcriptomic",
        label = sprintf("%s %s", chem, direction),
        source = normalizePath(path, winslash = "/", mustWork = TRUE),
        genes = normalize_ids(readr::read_lines(path))
      )
    }
  }

  literature_dir <- file.path(OUT_DIR, "literature_tf_expansion", "Literature_DEGs")
  if (dir.exists(literature_dir)) {
    paths <- sort(list.files(
      literature_dir,
      pattern = "(?:EXPANDED|ORIGINAL)\\.csv$",
      full.names = TRUE
    ))
    for (path in paths) {
      stem <- tools::file_path_sans_ext(basename(path))
      display_label <- sub("^TFEXP_", "", stem)
      display_label <- sub("_ALL", "", display_label, fixed = TRUE)
      display_label <- sub("_(EXPANDED|ORIGINAL)$", "", display_label)
      display_label <- gsub("_", " ", display_label, fixed = TRUE)
      key <- paste("Literature", stem, sep = "::")
      records[[key]] <- list(
        branch = "Literature",
        label = display_label,
        source = normalizePath(path, winslash = "/", mustWork = TRUE),
        genes = read_gene_column(path)
      )
    }
  }

  records
}

gene_sets <- discover_gene_sets()
if (!length(gene_sets)) {
  stop("No gene sets were discovered. Run steps 01 and 02 first.", call. = FALSE)
}

data("human_ens_clusters", package = "AOPfingerprintR")
data("human_ens_aop", package = "AOPfingerprintR")
families <- list(KE = human_ens_clusters, AOP = human_ens_aop)

prepare_annotation <- function(annotation_sets, background) {
  background <- normalize_ids(background)
  background_index <- setNames(seq_along(background), background)
  filtered <- lapply(annotation_sets, function(x) intersect(normalize_ids(x), background))
  keep <- lengths(filtered) > 0L & lengths(filtered) < length(background)
  filtered <- filtered[keep]
  if (!length(filtered)) stop("No annotation terms overlap the configured background.", call. = FALSE)

  term_index <- rep(seq_along(filtered), lengths(filtered))
  gene_index <- unname(background_index[unlist(filtered, use.names = FALSE)])
  membership <- Matrix::sparseMatrix(
    i = gene_index,
    j = term_index,
    x = 1,
    dims = c(length(background), length(filtered))
  )
  list(
    background = background,
    background_index = background_index,
    membership = membership,
    term_ids = names(filtered),
    term_sizes = lengths(filtered)
  )
}

annotation_cache <- lapply(names(backgrounds), function(branch) {
  lapply(families, prepare_annotation, background = backgrounds[[branch]])
})
names(annotation_cache) <- names(backgrounds)

ora_matrix <- function(overlap, term_sizes, universe_n, set_size) {
  iterations <- nrow(overlap)
  term_count <- ncol(overlap)
  pvalues <- stats::phyper(
    q = as.vector(overlap) - 1,
    m = rep(term_sizes, each = iterations),
    n = rep(universe_n - term_sizes, each = iterations),
    k = set_size,
    lower.tail = FALSE
  )
  pvalues <- matrix(pvalues, nrow = iterations, ncol = term_count)
  t(apply(pvalues, 1, stats::p.adjust, method = "BH"))
}

evaluate_observed <- function(genes, prepared, family, set_id, label, branch, source) {
  genes <- intersect(normalize_ids(genes), prepared$background)
  if (!length(genes)) return(NULL)
  indices <- unname(prepared$background_index[genes])
  overlap <- matrix(Matrix::colSums(prepared$membership[indices, , drop = FALSE]), nrow = 1)
  fdr <- ora_matrix(overlap, prepared$term_sizes, length(prepared$background), length(genes))
  pvalue <- stats::phyper(
    overlap[1, ] - 1,
    prepared$term_sizes,
    length(prepared$background) - prepared$term_sizes,
    length(genes),
    lower.tail = FALSE
  )
  term_table <- data.frame(
    set_id = set_id,
    label = label,
    branch = branch,
    family = family,
    term_id = prepared$term_ids,
    universe_size = length(prepared$background),
    gene_set_size = length(genes),
    term_size = prepared$term_sizes,
    overlap = as.integer(overlap[1, ]),
    pvalue = pvalue,
    fdr = as.numeric(fdr[1, ]),
    source = source,
    stringsAsFactors = FALSE
  )
  list(
    genes = genes,
    terms = term_table,
    metrics = data.frame(
      set_id = set_id,
      label = label,
      branch = branch,
      family = family,
      universe_size = length(prepared$background),
      gene_set_size = length(genes),
      n_significant_terms = sum(fdr[1, ] < FDR_CUTOFF, na.rm = TRUE),
      top_neglog10_fdr = -log10(max(min(fdr[1, ], na.rm = TRUE), .Machine$double.xmin)),
      source = source,
      stringsAsFactors = FALSE
    )
  )
}

simulate_null <- function(set_size, prepared, iterations, seed) {
  universe_n <- length(prepared$background)
  if (set_size > universe_n) stop("Gene-set size exceeds the background.", call. = FALSE)
  set.seed(seed)
  batches <- split(seq_len(iterations), ceiling(seq_len(iterations) / NULL_BATCH_SIZE))
  results <- lapply(batches, function(batch_ids) {
    batch_n <- length(batch_ids)
    sampled <- replicate(batch_n, sample.int(universe_n, set_size, replace = FALSE))
    draws <- Matrix::sparseMatrix(
      i = rep(seq_len(batch_n), each = set_size),
      j = as.vector(sampled),
      x = 1,
      dims = c(batch_n, universe_n)
    )
    overlap <- as.matrix(draws %*% prepared$membership)
    fdr <- ora_matrix(overlap, prepared$term_sizes, universe_n, set_size)
    data.frame(
      iteration = batch_ids,
      n_significant_terms = rowSums(fdr < FDR_CUTOFF, na.rm = TRUE),
      top_neglog10_fdr = -log10(pmax(apply(fdr, 1, min, na.rm = TRUE), .Machine$double.xmin))
    )
  })
  dplyr::bind_rows(results)
}

observed_rows <- list()
term_rows <- list()
null_rows <- list()
counter <- 0L

for (set_id in names(gene_sets)) {
  item <- gene_sets[[set_id]]
  for (family in names(families)) {
    counter <- counter + 1L
    prepared <- annotation_cache[[item$branch]][[family]]
    observed <- evaluate_observed(
      item$genes, prepared, family, set_id, item$label, item$branch, item$source
    )
    if (is.null(observed)) next
    message(sprintf(
      "[%d] %s | %s | n=%d",
      counter, item$label, family, observed$metrics$gene_set_size
    ))
    null <- simulate_null(
      observed$metrics$gene_set_size,
      prepared,
      NULL_ITERATIONS,
      NULL_SEED + counter * 1009L
    )
    null$set_id <- set_id
    null$label <- item$label
    null$branch <- item$branch
    null$family <- family

    observed_rows[[length(observed_rows) + 1L]] <- observed$metrics
    term_rows[[length(term_rows) + 1L]] <- observed$terms
    null_rows[[length(null_rows) + 1L]] <- null
  }
}

observed_df <- dplyr::bind_rows(observed_rows)
term_df <- dplyr::bind_rows(term_rows)
null_df <- dplyr::bind_rows(null_rows)

summary_df <- observed_df |>
  dplyr::left_join(
    null_df |>
      dplyr::group_by(set_id, family) |>
      dplyr::summarise(
        null_median = stats::median(n_significant_terms),
        null_q025 = stats::quantile(n_significant_terms, 0.025),
        null_q975 = stats::quantile(n_significant_terms, 0.975),
        null_mean = mean(n_significant_terms),
        null_sd = stats::sd(n_significant_terms),
        .groups = "drop"
      ),
    by = c("set_id", "family")
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    empirical_p = (1 + sum(
      null_df$set_id == set_id &
        null_df$family == family &
        null_df$n_significant_terms >= n_significant_terms
    )) / (NULL_ITERATIONS + 1),
    z_score = ifelse(is.finite(null_sd) && null_sd > 0,
                     (n_significant_terms - null_mean) / null_sd, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(family) |>
  dplyr::mutate(empirical_fdr = stats::p.adjust(empirical_p, method = "BH")) |>
  dplyr::ungroup()

readr::write_csv(observed_df, file.path(SPEC_DIR, "observed_enrichment_metrics.csv"))
readr::write_csv(term_df, file.path(SPEC_DIR, "observed_term_level_enrichment.csv"))
readr::write_csv(null_df, file.path(SPEC_DIR, "null_enrichment_metrics.csv.gz"))
readr::write_csv(summary_df, file.path(SPEC_DIR, "specificity_summary.csv"))

plot_df <- summary_df |>
  dplyr::mutate(
    display_label = ifelse(branch == "Transcriptomic", label, paste0("Lit: ", label)),
    display_label = factor(display_label, levels = rev(unique(display_label))),
    specificity = ifelse(empirical_fdr < 0.05, "Empirical FDR < 0.05", "Not significant")
  )

p <- ggplot(plot_df, aes(y = display_label)) +
  geom_errorbar(
    aes(xmin = null_q025, xmax = null_q975),
    orientation = "y", width = 0.18, linewidth = 0.7, colour = "#9AA3A8"
  ) +
  geom_point(aes(x = null_median), shape = 124, size = 5, colour = "#59636A") +
  geom_point(
    aes(x = n_significant_terms, fill = specificity),
    shape = 21, size = 3.4, stroke = 0.35, colour = "#1B2730"
  ) +
  facet_grid(branch ~ family, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Empirical FDR < 0.05" = "#D95F02", "Not significant" = "#F3C969")) +
  labs(
    x = sprintf("Significant terms (BH FDR < %.2f)", FDR_CUTOFF),
    y = NULL,
    fill = NULL,
    title = "Observed KE/AOP enrichment versus size-matched random gene sets",
    subtitle = sprintf("Point: observed; vertical tick and line: null median and 95%% interval (%s iterations)",
                       scales::comma(NULL_ITERATIONS))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

ggsave(file.path(SPEC_DIR, "Figure_AOP_specificity_null.png"), p,
       width = 12.5, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(SPEC_DIR, "Figure_AOP_specificity_null.pdf"), p,
       width = 12.5, height = 8.5, device = cairo_pdf, bg = "white")

provenance <- c(
  "KE/AOP specificity analysis",
  sprintf("Generated: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("Iterations: %d", NULL_ITERATIONS),
  sprintf("Seed: %d", NULL_SEED),
  sprintf("Directions: %s", paste(NULL_DIRECTIONS, collapse = ",")),
  sprintf("Term-level test: one-sided hypergeometric over-representation"),
  sprintf("Term-level family: each gene set x annotation family (KE or AOP)"),
  sprintf("Term-level correction: Benjamini-Hochberg; FDR cutoff %.3f", FDR_CUTOFF),
  "Empirical statistic: number of significant terms",
  "Empirical p: (1 + null values >= observed) / (iterations + 1)",
  "Empirical correction: Benjamini-Hochberg across evaluated sets within KE and AOP families",
  sprintf("Transcriptomic background: %s", normalizePath(BACKGROUND_GSE_FILE, winslash = "/")),
  sprintf("Literature background: %s", normalizePath(BACKGROUND_LIT_FILE, winslash = "/"))
)
writeLines(provenance, file.path(SPEC_DIR, "specificity_provenance.txt"))

cat("Specificity analysis complete:", normalizePath(SPEC_DIR, winslash = "/"), "\n")

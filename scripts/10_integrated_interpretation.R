#!/usr/bin/env Rscript
# Integrated interpretation of functional enrichment, AOP fingerprints,
# transcription-factor expansion, and optional literature-signature projection.
# All inputs are discovered from standard pipeline outputs.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "readr", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap", "scales",
  "stringr", "forcats"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

RESULT_DIR <- file.path(OUT_DIR, "Interpretation", "Integrated")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N <- suppressWarnings(as.integer(
  Sys.getenv("AOP_PUBLICATION_TOP_N", unset = "10")
))
if (!is.finite(TOP_N) || TOP_N < 1L) TOP_N <- 10L
TOP_AOPS <- suppressWarnings(as.integer(
  Sys.getenv("AOP_PUBLICATION_TOP_AOPS", unset = "25")
))
if (!is.finite(TOP_AOPS) || TOP_AOPS < 1L) TOP_AOPS <- 25L

PUBLICATION_DPI <- suppressWarnings(as.integer(
  Sys.getenv("AOP_PUBLICATION_DPI", unset = "600")
))
if (!is.finite(PUBLICATION_DPI) || PUBLICATION_DPI < 72L) PUBLICATION_DPI <- 600L

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

read_csv_quiet <- function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
}

first_nonmissing <- function(x, fallback = "Uncategorized") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) sort(unique(x))[[1]] else fallback
}

classify_source <- function(experiment) {
  ifelse(
    grepl("^TFEXP_", experiment, ignore.case = TRUE),
    "Literature + TF expansion",
    "Transcriptomic"
  )
}

parse_tf_experiment <- function(experiment) {
  match <- stringr::str_match(
    experiment,
    "^TFEXP_(.+?)_ALL_(AB|ABC|ORIG)_(EXPANDED|ORIGINAL)$"
  )
  tibble::tibble(
    Experiment = experiment,
    Chemical = match[, 2],
    Level = dplyr::recode(match[, 3], ORIG = "Original", AB = "AB", ABC = "ABC")
  )
}

strip_plot_prefix <- function(x) {
  stringr::str_wrap(stringr::str_remove(x, "^.*___"), width = 48)
}

pretty_experiment <- function(experiment) {
  result <- as.character(experiment)
  matched <- stringr::str_match(
    result,
    "^TFEXP_(.+?)_ALL_(AB|ABC|ORIG)_(EXPANDED|ORIGINAL)$"
  )
  is_tf <- !is.na(matched[, 2])
  level <- c(ORIG = "Original", AB = "AB", ABC = "ABC")[matched[, 3]]
  result[is_tf] <- paste(matched[is_tf, 2], level[is_tf], sep = " · ")
  result
}

level_order <- c("Original", "AB", "ABC")
level_colors <- c(Original = "#666666", AB = "#1B9E77", ABC = "#D95F02")
source_colors <- c(
  "Transcriptomic" = "#3B6FB6",
  "Literature + TF expansion" = "#D55E00"
)

report_lines <- c(
  "# Integrated interpretation summary",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  ""
)

# =============================================================================
# Functional-enrichment recurrence and landscape
# =============================================================================
summary_dir <- file.path(ENR_DIR, "_summaries")
functional_files <- list.files(
  summary_dir,
  pattern = "^merged_.*_significant\\.csv$",
  full.names = TRUE
)

functional <- NULL
functional_recurrence <- NULL
if (length(functional_files)) {
  functional <- bind_rows(lapply(functional_files, function(path) {
    data <- read_csv_quiet(path)
    required <- c("Description", "Experiment", "p.adjust")
    missing <- setdiff(required, names(data))
    if (length(missing)) {
      warning(basename(path), " skipped; missing: ", paste(missing, collapse = ", "))
      return(NULL)
    }
    data$Collection <- sub(
      "^merged_(.*)_significant\\.csv$",
      "\\1",
      basename(path)
    )
    data
  }))

  if (nrow(functional)) {
    if (!"Count" %in% names(functional)) functional$Count <- NA_real_
    functional <- functional %>%
      mutate(
        Description = as.character(Description),
        Experiment = as.character(Experiment),
        Collection = as.character(Collection),
        p.adjust = as.numeric(p.adjust),
        Count = as.numeric(Count),
        neglog10_padj = -log10(pmax(p.adjust, .Machine$double.xmin)),
        Source = classify_source(Experiment)
      ) %>%
      filter(!is.na(Description), nzchar(Description), is.finite(neglog10_padj))

    functional_recurrence <- functional %>%
      group_by(Collection, Description) %>%
      summarise(
        n_experiments = n_distinct(Experiment),
        n_transcriptomic = n_distinct(Experiment[Source == "Transcriptomic"]),
        n_literature_tf = n_distinct(Experiment[Source == "Literature + TF expansion"]),
        best_padj = min(p.adjust, na.rm = TRUE),
        max_neglog10_padj = max(neglog10_padj, na.rm = TRUE),
        median_neglog10_padj = median(neglog10_padj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(Collection, desc(n_experiments), desc(max_neglog10_padj))
    readr::write_csv(
      functional_recurrence,
      file.path(RESULT_DIR, "functional_term_recurrence.csv")
    )

    top_functional <- functional_recurrence %>%
      group_by(Collection) %>%
      slice_max(
        order_by = n_experiments * 1e6 + pmin(max_neglog10_padj, 1e5),
        n = TOP_N,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      arrange(Collection, n_experiments, max_neglog10_padj) %>%
      mutate(
        Plot_term = factor(
          paste(Collection, Description, sep = "___"),
          levels = unique(paste(Collection, Description, sep = "___"))
        )
      )

    p_functional_recurrence <- ggplot(
      top_functional,
      aes(n_experiments, Plot_term, color = max_neglog10_padj)
    ) +
      geom_segment(aes(x = 0, xend = n_experiments, yend = Plot_term), color = "grey80") +
      geom_point(size = 3.2) +
      facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
      scale_y_discrete(labels = strip_plot_prefix) +
      scale_color_viridis_c(option = "C", trans = "sqrt") +
      scale_x_continuous(breaks = scales::breaks_pretty()) +
      labs(
        title = "Recurrent functional terms across experiments",
        subtitle = sprintf("Top %d terms per collection", TOP_N),
        x = "Experiments with significant enrichment",
        y = NULL,
        color = expression(max(-log[10](FDR)))
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    save_plot(
      p_functional_recurrence,
      "Figure_functional_term_recurrence",
      width = 11,
      height = 16
    )

    selected_terms <- top_functional %>% select(Collection, Description)
    functional_landscape <- functional %>%
      inner_join(selected_terms, by = c("Collection", "Description")) %>%
      group_by(Collection, Description, Experiment, Source) %>%
      summarise(
        neglog10_padj = max(neglog10_padj, na.rm = TRUE),
        Count = max(Count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Count = ifelse(is.finite(Count), Count, 1),
        Experiment = factor(
          Experiment,
          levels = unique(c(
            sort(Experiment[Source == "Transcriptomic"]),
            sort(Experiment[Source == "Literature + TF expansion"])
          ))
        ),
        Plot_term = paste(Collection, Description, sep = "___")
      )
    readr::write_csv(
      functional_landscape,
      file.path(RESULT_DIR, "functional_enrichment_landscape.csv")
    )

    p_functional_landscape <- ggplot(
      functional_landscape,
      aes(Experiment, Plot_term, size = Count, color = neglog10_padj)
    ) +
      geom_point(alpha = 0.85) +
      facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
      scale_y_discrete(labels = strip_plot_prefix) +
      scale_size_continuous(range = c(1.2, 7), trans = "sqrt") +
      scale_color_viridis_c(option = "B", trans = "sqrt") +
      labs(
        title = "Functional-enrichment landscape",
        subtitle = "Bubble size: enriched genes; color: enrichment significance",
        x = NULL,
        y = NULL,
        size = "Gene count",
        color = expression(-log[10](FDR))
      ) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    save_plot(
      p_functional_landscape,
      "Figure_functional_enrichment_landscape",
      width = 18,
      height = 18
    )

    top_by_collection <- functional_recurrence %>%
      group_by(Collection) %>%
      slice_head(n = 1) %>%
      ungroup()
    report_lines <- c(
      report_lines,
      "## Functional enrichment",
      "",
      sprintf("- Significant merged terms analyzed: %s.", scales::comma(nrow(functional))),
      vapply(seq_len(nrow(top_by_collection)), function(i) {
        sprintf(
          "- %s: most recurrent term is `%s` (%d experiments).",
          top_by_collection$Collection[[i]],
          top_by_collection$Description[[i]],
          top_by_collection$n_experiments[[i]]
        )
      }, character(1)),
      ""
    )
  }
} else {
  message("Functional summaries absent; skipping functional synthesis.")
}

# =============================================================================
# AOP recurrence, annotations, and experiment similarity
# =============================================================================
aop_files <- c(
  GSE = file.path(AOP_DIR, "AOP_fingerprint_GSE_enriched.csv"),
  Literature = file.path(AOP_DIR, "AOP_fingerprint_Literature_enriched.csv")
)
aop_files <- aop_files[file.exists(aop_files)]

aop_scores <- NULL
aop_recurrence <- NULL
if (length(aop_files)) {
  aop_raw <- bind_rows(lapply(names(aop_files), function(input_source) {
    data <- read_csv_quiet(aop_files[[input_source]])
    data$Input_source <- input_source
    data
  }))
  required <- c("Experiment", "AOP_name", "neglog10_padj")
  missing <- setdiff(required, names(aop_raw))
  if (length(missing)) {
    warning("AOP synthesis skipped; missing columns: ", paste(missing, collapse = ", "))
  } else {
    for (column in c("SSbD_category", "Organ", "Endpoint", "Ke", "Ke_description")) {
      if (!column %in% names(aop_raw)) aop_raw[[column]] <- NA_character_
    }

    aop_scores <- aop_raw %>%
      mutate(
        Experiment = as.character(Experiment),
        AOP_name = as.character(AOP_name),
        neglog10_padj = as.numeric(neglog10_padj),
        Source = classify_source(Experiment)
      ) %>%
      filter(!is.na(AOP_name), nzchar(AOP_name), is.finite(neglog10_padj)) %>%
      group_by(Experiment, Source, AOP_name) %>%
      summarise(
        neglog10_padj = max(neglog10_padj, na.rm = TRUE),
        SSbD_category = first_nonmissing(SSbD_category),
        Organ = first_nonmissing(Organ),
        Endpoint = first_nonmissing(Endpoint),
        n_key_events = n_distinct(Ke, na.rm = TRUE),
        .groups = "drop"
      )

    aop_recurrence <- aop_scores %>%
      group_by(AOP_name) %>%
      summarise(
        n_experiments = n_distinct(Experiment),
        n_sources = n_distinct(Source),
        n_transcriptomic = n_distinct(Experiment[Source == "Transcriptomic"]),
        n_literature_tf = n_distinct(Experiment[Source == "Literature + TF expansion"]),
        max_neglog10_padj = max(neglog10_padj, na.rm = TRUE),
        median_neglog10_padj = median(neglog10_padj, na.rm = TRUE),
        SSbD_category = first_nonmissing(SSbD_category),
        Organ = first_nonmissing(Organ),
        Endpoint = first_nonmissing(Endpoint),
        .groups = "drop"
      ) %>%
      arrange(desc(n_experiments), desc(n_sources), desc(max_neglog10_padj))
    readr::write_csv(aop_recurrence, file.path(RESULT_DIR, "AOP_recurrence.csv"))

    top_aops <- aop_recurrence %>% slice_head(n = min(TOP_AOPS, nrow(aop_recurrence)))
    aop_landscape <- aop_scores %>%
      filter(AOP_name %in% top_aops$AOP_name) %>%
      mutate(
        AOP_name = factor(AOP_name, levels = rev(top_aops$AOP_name)),
        Experiment = factor(
          Experiment,
          levels = unique(c(
            sort(Experiment[Source == "Transcriptomic"]),
            sort(Experiment[Source == "Literature + TF expansion"])
          ))
        )
      )
    readr::write_csv(aop_landscape, file.path(RESULT_DIR, "AOP_landscape.csv"))

    p_aop <- ggplot(
      aop_landscape,
      aes(Experiment, AOP_name, size = n_key_events, color = neglog10_padj)
    ) +
      geom_point(alpha = 0.85) +
      scale_size_continuous(range = c(1.5, 7)) +
      scale_color_viridis_c(option = "C", trans = "sqrt") +
      labs(
        title = "Recurrent Adverse Outcome Pathway landscape",
        subtitle = "Top recurrent AOPs across transcriptomic and literature-derived experiments",
        x = NULL,
        y = NULL,
        size = "Enriched KEs",
        color = expression(-log[10](FDR))
      ) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    aop_height <- max(6, min(20, 3 + 0.36 * nrow(top_aops)))
    save_plot(p_aop, "Figure_AOP_landscape", width = 18, height = aop_height)

    annotation_long <- aop_scores %>%
      select(Experiment, Source, AOP_name, SSbD_category, Organ, Endpoint) %>%
      pivot_longer(
        c(SSbD_category, Organ, Endpoint),
        names_to = "Dimension",
        values_to = "Category"
      ) %>%
      filter(!is.na(Category), nzchar(Category), Category != "Uncategorized") %>%
      distinct(Source, Experiment, AOP_name, Dimension, Category) %>%
      count(Dimension, Source, Category, name = "n_AOP_experiment_pairs") %>%
      group_by(Dimension) %>%
      mutate(total = sum(n_AOP_experiment_pairs)) %>%
      ungroup()

    top_annotations <- annotation_long %>%
      group_by(Dimension, Category) %>%
      summarise(total = sum(n_AOP_experiment_pairs), .groups = "drop") %>%
      group_by(Dimension) %>%
      slice_max(total, n = TOP_N, with_ties = FALSE) %>%
      ungroup() %>%
      select(Dimension, Category)
    annotation_plot_data <- annotation_long %>%
      inner_join(top_annotations, by = c("Dimension", "Category")) %>%
      mutate(Plot_category = paste(Dimension, Category, sep = "___"))
    readr::write_csv(
      annotation_long,
      file.path(RESULT_DIR, "AOP_annotation_counts.csv")
    )

    p_annotations <- ggplot(
      annotation_plot_data,
      aes(n_AOP_experiment_pairs, forcats::fct_reorder(Plot_category, n_AOP_experiment_pairs), fill = Source)
    ) +
      geom_col(position = "dodge") +
      facet_grid(Dimension ~ ., scales = "free_y", space = "free_y") +
      scale_y_discrete(labels = strip_plot_prefix) +
      scale_fill_manual(values = source_colors) +
      labs(
        title = "AOP interpretation by hazard category, organ, and endpoint",
        x = "AOP-experiment pairs",
        y = NULL,
        fill = "Evidence source"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "top"
      )
    save_plot(p_annotations, "Figure_AOP_annotation_profile", width = 11, height = 12)

    aop_wide <- aop_scores %>%
      select(AOP_name, Experiment, neglog10_padj) %>%
      pivot_wider(names_from = Experiment, values_from = neglog10_padj, values_fill = 0)
    aop_matrix <- as.matrix(aop_wide[, -1, drop = FALSE])
    if (ncol(aop_matrix) > 1L) {
      similarity <- stats::cor(aop_matrix, method = "spearman", use = "pairwise.complete.obs")
      similarity[is.na(similarity)] <- 0
      diag(similarity) <- 1
      readr::write_csv(
        as.data.frame(similarity) %>% tibble::rownames_to_column("Experiment"),
        file.path(RESULT_DIR, "AOP_experiment_similarity.csv")
      )
      plot_similarity <- similarity
      dimnames(plot_similarity) <- list(
        pretty_experiment(rownames(plot_similarity)),
        pretty_experiment(colnames(plot_similarity))
      )
      pheatmap::pheatmap(
        plot_similarity,
        color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
        breaks = seq(-1, 1, length.out = 102),
        border_color = "white",
        main = "Experiment similarity based on AOP fingerprints",
        filename = file.path(RESULT_DIR, "Figure_AOP_similarity.png"),
        width = 9,
        height = 8
      )
      pheatmap::pheatmap(
        plot_similarity,
        color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
        breaks = seq(-1, 1, length.out = 102),
        border_color = "white",
        main = "Experiment similarity based on AOP fingerprints",
        filename = file.path(RESULT_DIR, "Figure_AOP_similarity.pdf"),
        width = 9,
        height = 8
      )
    }

    report_lines <- c(
      report_lines,
      "## Adverse Outcome Pathways",
      "",
      sprintf("- AOP-experiment associations analyzed: %s.", scales::comma(nrow(aop_scores))),
      sprintf(
        "- Most recurrent AOP: `%s` (%d experiments; category: %s; organ: %s).",
        aop_recurrence$AOP_name[[1]],
        aop_recurrence$n_experiments[[1]],
        aop_recurrence$SSbD_category[[1]],
        aop_recurrence$Organ[[1]]
      ),
      ""
    )
  }
} else {
  message("AOP fingerprint tables absent; skipping AOP synthesis.")
}

# =============================================================================
# TF-expansion impact on gene sets and downstream yield
# =============================================================================
tf_root <- file.path(OUT_DIR, "literature_tf_expansion")
tf_gene_files <- if (dir.exists(tf_root)) {
  list.files(
    tf_root,
    pattern = "^TFEXP_.*_(AB|ABC|ORIG)_(EXPANDED|ORIGINAL)\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
} else {
  character()
}

tf_sizes <- NULL
if (length(tf_gene_files)) {
  tf_sizes <- bind_rows(lapply(tf_gene_files, function(path) {
    parsed <- stringr::str_match(
      tools::file_path_sans_ext(basename(path)),
      "^TFEXP_(.+?)_ALL_(AB|ABC|ORIG)_(EXPANDED|ORIGINAL)$"
    )
    if (is.na(parsed[1, 2])) return(NULL)
    data <- read_csv_quiet(path)
    gene_col <- intersect(c("Gene", "Feature", "ENSEMBL", "SYMBOL"), names(data))
    if (!length(gene_col)) return(NULL)
    tibble::tibble(
      Chemical = parsed[1, 2],
      Level = dplyr::recode(parsed[1, 3], ORIG = "Original", AB = "AB", ABC = "ABC"),
      n_genes = n_distinct(data[[gene_col[[1]]]])
    )
  })) %>%
    distinct(Chemical, Level, .keep_all = TRUE) %>%
    mutate(Level = factor(Level, levels = level_order)) %>%
    arrange(Chemical, Level)
  readr::write_csv(tf_sizes, file.path(RESULT_DIR, "TF_expansion_gene_set_sizes.csv"))

  p_tf_sizes <- ggplot(tf_sizes, aes(Chemical, n_genes, fill = Level)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(
      aes(label = scales::comma(n_genes)),
      position = position_dodge(width = 0.8),
      vjust = -0.25,
      size = 3
    ) +
    scale_fill_manual(values = level_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Impact of transcription-factor expansion on gene-set size",
      x = NULL,
      y = "Unique genes",
      fill = "Gene set"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  save_plot(p_tf_sizes, "Figure_TF_expansion_gene_set_size", width = 8, height = 5)

  tf_file_chemicals <- stringr::str_match(
    tools::file_path_sans_ext(basename(tf_gene_files)),
    "^TFEXP_(.+?)_ALL_"
  )[, 2]
  tf_jaccard <- bind_rows(lapply(split(tf_gene_files, tf_file_chemicals), function(paths) {
    parsed <- lapply(paths, function(path) {
      name_match <- stringr::str_match(
        tools::file_path_sans_ext(basename(path)),
        "^TFEXP_(.+?)_ALL_(AB|ABC|ORIG)_(EXPANDED|ORIGINAL)$"
      )
      data <- read_csv_quiet(path)
      gene_col <- intersect(c("Gene", "Feature", "ENSEMBL", "SYMBOL"), names(data))
      if (!length(gene_col)) return(NULL)
      gene_col <- gene_col[[1]]
      list(
        Chemical = name_match[1, 2],
        Level = dplyr::recode(name_match[1, 3], ORIG = "Original", AB = "AB", ABC = "ABC"),
        Genes = unique(na.omit(as.character(data[[gene_col]])))
      )
    })
    parsed <- Filter(Negate(is.null), parsed)
    if (length(parsed) < 2L) return(NULL)
    pairs <- utils::combn(seq_along(parsed), 2, simplify = FALSE)
    bind_rows(lapply(pairs, function(pair) {
      a <- parsed[[pair[[1]]]]
      b <- parsed[[pair[[2]]]]
      union_size <- length(union(a$Genes, b$Genes))
      tibble::tibble(
        Chemical = a$Chemical,
        Level_A = a$Level,
        Level_B = b$Level,
        n_A = length(a$Genes),
        n_B = length(b$Genes),
        n_intersection = length(intersect(a$Genes, b$Genes)),
        Jaccard = ifelse(union_size > 0, length(intersect(a$Genes, b$Genes)) / union_size, NA_real_)
      )
    }))
  }))
  readr::write_csv(tf_jaccard, file.path(RESULT_DIR, "TF_expansion_gene_set_similarity.csv"))

  if (!is.null(functional) && nrow(functional)) {
    tf_map <- parse_tf_experiment(unique(functional$Experiment))
    tf_functional_yield <- functional %>%
      inner_join(tf_map, by = "Experiment") %>%
      filter(!is.na(Chemical), !is.na(Level)) %>%
      group_by(Chemical, Level, Collection) %>%
      summarise(n_significant_terms = n_distinct(Description), .groups = "drop") %>%
      mutate(Level = factor(Level, levels = level_order))
    readr::write_csv(
      tf_functional_yield,
      file.path(RESULT_DIR, "TF_expansion_functional_yield.csv")
    )

    p_tf_functional <- ggplot(
      tf_functional_yield,
      aes(Collection, n_significant_terms, fill = Level)
    ) +
      geom_col(position = "dodge") +
      facet_wrap(~Chemical, scales = "free_y") +
      scale_fill_manual(values = level_colors, drop = FALSE) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = "Effect of TF expansion on functional-enrichment yield",
        x = NULL,
        y = "Significant terms",
        fill = "Gene set"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "top"
      )
    save_plot(p_tf_functional, "Figure_TF_expansion_functional_yield", width = 11, height = 5.5)
  }

  if (!is.null(aop_scores) && nrow(aop_scores)) {
    tf_aop_map <- parse_tf_experiment(unique(aop_scores$Experiment))
    tf_aop_yield <- aop_scores %>%
      inner_join(tf_aop_map, by = "Experiment") %>%
      filter(!is.na(Chemical), !is.na(Level)) %>%
      group_by(Chemical, Level) %>%
      summarise(
        n_AOPs = n_distinct(AOP_name),
        n_AOP_experiment_pairs = n(),
        median_neglog10_padj = median(neglog10_padj),
        .groups = "drop"
      ) %>%
      mutate(Level = factor(Level, levels = level_order))
    readr::write_csv(tf_aop_yield, file.path(RESULT_DIR, "TF_expansion_AOP_yield.csv"))

    p_tf_aop <- ggplot(tf_aop_yield, aes(Chemical, n_AOPs, fill = Level)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = level_colors, drop = FALSE) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = "Effect of TF expansion on AOP yield",
        x = NULL,
        y = "Distinct enriched AOPs",
        fill = "Gene set"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")
    save_plot(p_tf_aop, "Figure_TF_expansion_AOP_yield", width = 8, height = 5)
  }

  tf_fold <- tf_sizes %>%
    mutate(Level = as.character(Level)) %>%
    pivot_wider(names_from = Level, values_from = n_genes)
  for (column in level_order) {
    if (!column %in% names(tf_fold)) tf_fold[[column]] <- NA_real_
  }
  tf_fold <- tf_fold %>%
    mutate(
      AB_fold = ifelse(is.finite(Original) & Original > 0, AB / Original, NA_real_),
      ABC_fold = ifelse(is.finite(Original) & Original > 0, ABC / Original, NA_real_)
    )
  report_lines <- c(
    report_lines,
    "## Transcription-factor expansion",
    "",
    vapply(seq_len(nrow(tf_fold)), function(i) {
      sprintf(
        "- %s: AB expands the original set %.2fx; ABC expands it %.2fx.",
        tf_fold$Chemical[[i]], tf_fold$AB_fold[[i]], tf_fold$ABC_fold[[i]]
      )
    }, character(1)),
    ""
  )
} else {
  message("TF-expansion tables absent; skipping TF-impact synthesis.")
}

# =============================================================================
# Literature-signature projection dot plot (optional step 08)
# =============================================================================
projection_file <- file.path(
  OUT_DIR,
  "Signature_Projection",
  "fgsea_all_contrasts_literature_sets.csv"
)
if (file.exists(projection_file)) {
  projection <- read_csv_quiet(projection_file) %>%
    mutate(
      NES = as.numeric(NES),
      padj = as.numeric(padj),
      neglog10_padj = -log10(pmax(padj, .Machine$double.xmin))
    ) %>%
    filter(is.finite(NES), is.finite(neglog10_padj))
  readr::write_csv(projection, file.path(RESULT_DIR, "signature_projection_plot_data.csv"))

  if (nrow(projection)) {
    limit <- max(abs(projection$NES), na.rm = TRUE)
    p_projection <- ggplot(
      projection,
      aes(Contrast, pathway, color = NES, size = neglog10_padj)
    ) +
      geom_point(alpha = 0.9) +
      scale_color_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-limit, limit)
      ) +
      scale_size_continuous(range = c(2, 8)) +
      labs(
        title = "Literature signatures projected onto transcriptomic contrasts",
        x = NULL,
        y = "Literature signature",
        color = "NES",
        size = expression(-log[10](FDR))
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold")
      )
    save_plot(p_projection, "Figure_signature_projection", width = 8.5, height = 5)

    significant_projection <- projection %>% filter(padj < FDR_CUTOFF)
    report_lines <- c(
      report_lines,
      "## Literature-signature projection",
      "",
      sprintf(
        "- Significant signature-contrast associations at FDR < %.2g: %d.",
        FDR_CUTOFF,
        nrow(significant_projection)
      ),
      ""
    )
  }
}

report_lines <- c(
  report_lines,
  "## Interpretation note",
  "",
  "These analyses rank recurrence, concordance, and enrichment strength to aid interpretation.",
  "They are exploratory summaries and should be evaluated with study design, dose, time, and biological context."
)
writeLines(report_lines, file.path(RESULT_DIR, "interpretation_summary.md"), useBytes = TRUE)

message("Integrated interpretation complete: ", normalizePath(RESULT_DIR, winslash = "/"))

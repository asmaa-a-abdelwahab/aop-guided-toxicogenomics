#!/usr/bin/env Rscript
# Assemble journal-size static figures from pipeline and app-export evidence.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "readr", "dplyr", "tidyr", "ggplot2", "scales", "stringr",
  "patchwork", "jsonlite", "igraph", "ggraph", "ggrepel", "tibble"
))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(igraph)
  library(ggraph)
})

PUB_DIR <- file.path(OUT_DIR, "Publication_figures")
dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)
PUB_DPI <- as.integer(numeric_setting("AOP_PUBLICATION_DPI", 600))

palette <- c(
  teal = "#0E7C7B", orange = "#E4572E", gold = "#F3A712",
  navy = "#1D3557", blue = "#457B9D", pale = "#EAF4F4",
  grey = "#687078", lightgrey = "#D8DEE3"
)

theme_publication <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, colour = "#17232B"),
      plot.subtitle = element_text(size = base_size, colour = "#46535C"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "#27333B"),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

manifest <- list()
save_publication_plot <- function(plot, stem, width, height, sources) {
  png_path <- file.path(PUB_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(PUB_DIR, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = PUB_DPI, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  manifest[[length(manifest) + 1L]] <<- data.frame(
    figure = stem,
    png = normalizePath(png_path, winslash = "/", mustWork = TRUE),
    pdf = normalizePath(pdf_path, winslash = "/", mustWork = TRUE),
    source_files = paste(normalizePath(sources, winslash = "/", mustWork = TRUE), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

must_read_csv <- function(path) {
  if (!file.exists(path)) stop("Required figure source not found: ", path, call. = FALSE)
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
}

deg_dir <- file.path(OUT_DIR, "Interpretation", "DEG_concordance")
integrated_dir <- file.path(OUT_DIR, "Interpretation", "Integrated")

# Figure 1: workflow ------------------------------------------------------------
workflow_nodes <- data.frame(
  x = c(1, 1, 3, 3, 5, 5, 7),
  y = c(2.6, 0.7, 2.6, 0.7, 2.6, 0.7, 1.65),
  label = c(
    "GEO FASTQ\nFastQC + MultiQC", "Curated literature\nsearch + signatures",
    "Counts + metadata\nvalidation / DESeq2", "DoRothEA\nA/B and A/B/C",
    "GO / KEGG / Reactome\nenrichment", "KE / AOP enrichment\nmatched-null specificity",
    "Mechanistic synthesis\nfigures + reusable outputs"
  ),
  branch = c("Experimental", "Literature", "Experimental", "Literature", "Shared", "Shared", "Output")
)
workflow_edges <- data.frame(
  x = c(1.65, 1.65, 3.65, 3.65, 5.65, 5.65),
  y = c(2.6, 0.7, 2.6, 0.7, 2.6, 0.7),
  xend = c(2.35, 2.35, 4.35, 4.35, 6.35, 6.35),
  yend = c(2.6, 0.7, 2.6, 0.7, 1.9, 1.4)
)
workflow_fill <- c(Experimental = "#D9F0EF", Literature = "#FDE6D8", Shared = "#E7EDF5", Output = "#FFF2C8")
p1 <- ggplot() +
  geom_segment(
    data = workflow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.9, colour = "#58666F",
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x = x, y = y, label = label, fill = branch),
    size = 3.7, fontface = "bold", linewidth = 0.35, label.padding = grid::unit(0.25, "lines"),
    colour = "#14212A", lineheight = 1.05
  ) +
  scale_fill_manual(values = workflow_fill) +
  annotate("text", x = 1, y = 3.35, label = "Measured branch", fontface = "bold", colour = palette[["teal"]], size = 4) +
  annotate("text", x = 1, y = -0.05, label = "Literature branch", fontface = "bold", colour = palette[["orange"]], size = 4) +
  coord_cartesian(xlim = c(0.1, 7.9), ylim = c(-0.3, 3.6), clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  labs(title = "AOP-guided toxicogenomics workflow") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), legend.position = "none")
save_publication_plot(p1, "Figure1_workflow", 13, 5.2, c(file.path(ROOT, "run_pipeline.R")))

# Figure 2: differential-expression burden ------------------------------------
burden_path <- file.path(deg_dir, "deg_burden_by_experiment.csv")
burden <- must_read_csv(burden_path) |>
  mutate(Experiment = factor(Experiment, levels = rev(intersect(CHEMS, Experiment)))) |>
  select(Experiment, n_up, n_down) |>
  pivot_longer(c(n_up, n_down), names_to = "Direction", values_to = "Genes") |>
  mutate(Direction = recode(Direction, n_up = "Upregulated", n_down = "Downregulated"))
p2 <- ggplot(burden, aes(x = Genes, y = Experiment, fill = Direction)) +
  geom_col(position = "stack", width = 0.68) +
  geom_text(
    data = burden |> group_by(Experiment) |> summarise(Genes = sum(Genes), .groups = "drop"),
    aes(x = Genes, y = Experiment, label = scales::comma(Genes)),
    hjust = -0.12, size = 3.6, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(Upregulated = palette[["orange"]], Downregulated = palette[["blue"]])) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Significant genes", y = NULL, fill = NULL,
    title = "Differential-expression burden across MXenes and positive controls",
    subtitle = sprintf("Benjamini-Hochberg FDR < %.2f and |log2 fold change| >= %.1f", FDR_CUTOFF, LFC_MIN)
  ) +
  theme_publication(11) +
  theme(panel.grid.major.y = element_blank())
save_publication_plot(p2, "Figure2_DEG_burden", 9.0, 5.2, burden_path)

# Figure 3: readable Jaccard heatmaps ------------------------------------------
overlap_path <- file.path(deg_dir, "pairwise_DEG_overlap.csv")
overlap <- must_read_csv(overlap_path) |>
  filter(Direction %in% c("ALL", "UP", "DOWN")) |>
  mutate(
    Experiment_A = factor(Experiment_A, levels = intersect(CHEMS, unique(Experiment_A))),
    Experiment_B = factor(Experiment_B, levels = rev(intersect(CHEMS, unique(Experiment_B)))),
    Direction = factor(Direction, levels = c("ALL", "UP", "DOWN")),
    label = ifelse(Experiment_A == rev(as.character(Experiment_B)), "", sprintf("%.2f", Jaccard))
  )
p3 <- ggplot(overlap, aes(x = Experiment_A, y = Experiment_B, fill = Jaccard)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3.5, fontface = "bold") +
  facet_wrap(~Direction, nrow = 1) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"), limits = c(0, 1)) +
  coord_fixed() +
  labs(
    x = NULL, y = NULL, fill = "Jaccard",
    title = "Pairwise overlap of significant gene sets",
    subtitle = "Separate panels retain direction while using a common 0-1 scale"
  ) +
  theme_publication(10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9.5),
    axis.text.y = element_text(size = 9.5),
    legend.position = "right"
  )
save_publication_plot(p3, "Figure3_DEG_overlap_readable", 13, 5.2, overlap_path)

# Figure 4: compact functional recurrence -------------------------------------
functional_path <- file.path(integrated_dir, "functional_term_recurrence.csv")
functional <- must_read_csv(functional_path) |>
  group_by(Collection) |>
  slice_max(order_by = n_experiments, n = 4, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    Description = stringr::str_wrap(Description, 30),
    term_key = paste(Collection, Description, sep = "::")
  ) |>
  arrange(Collection, n_experiments) |>
  mutate(term_key = factor(term_key, levels = unique(term_key)))
p4 <- ggplot(functional, aes(x = n_experiments, y = term_key, fill = max_neglog10_padj)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n_experiments), hjust = -0.25, size = 3.2, fontface = "bold") +
  facet_wrap(~Collection, scales = "free_y", ncol = 3, labeller = label_value) +
  scale_y_discrete(labels = function(x) sub("^[^:]+::", "", x)) +
  scale_x_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = c(0, 0.18))) +
  scale_fill_viridis_c(option = "C", end = 0.9) +
  labs(
    x = "Gene sets with significant enrichment", y = NULL,
    fill = "Maximum\n-log10(FDR)",
    title = "Recurrent functional programs across experimental and literature gene sets",
    subtitle = "Four most recurrent terms per ontology/pathway collection"
  ) +
  theme_publication(10) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8.7), legend.position = "right")
save_publication_plot(p4, "Figure4_functional_recurrence_readable", 13, 8.0, functional_path)

# Figure 5: TF expansion and signature projection -----------------------------
tf_path <- file.path(integrated_dir, "TF_expansion_gene_set_sizes.csv")
projection_path <- file.path(integrated_dir, "signature_projection_plot_data.csv")
tf_sizes <- must_read_csv(tf_path) |>
  mutate(Level = factor(Level, levels = c("Original", "AB", "ABC")))
tf_panel <- ggplot(tf_sizes, aes(x = Level, y = n_genes, fill = Chemical, group = Chemical)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = n_genes), position = position_dodge(width = 0.72), vjust = -0.35, size = 3.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  scale_fill_manual(values = c("HN-Ti3C2" = palette[["teal"]], "Ti3C2Tx" = palette[["orange"]])) +
  labs(x = NULL, y = "Unique genes", fill = NULL, title = "A  TF-regulon expansion") +
  theme_publication(10) + theme(panel.grid.major.x = element_blank())

projection <- must_read_csv(projection_path) |>
  mutate(
    pathway = stringr::str_replace_all(pathway, "_", " "),
    Contrast = factor(Contrast, levels = intersect(CHEMS, unique(Contrast))),
    significant = padj < FDR_CUTOFF
  )
projection_panel <- ggplot(projection, aes(x = Contrast, y = pathway)) +
  geom_point(aes(size = neglog10_padj, fill = NES, alpha = significant), shape = 21, colour = "#25323A", stroke = 0.3) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = NULL, y = NULL, fill = "NES", size = "-log10(FDR)", title = "B  Rank-based literature-signature projection") +
  theme_publication(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.minor = element_blank())
p5 <- tf_panel + projection_panel + plot_layout(widths = c(0.85, 1.35), guides = "collect") &
  theme(legend.position = "bottom")
save_publication_plot(p5, "Figure5_TF_expansion_signature_projection", 13, 5.7, c(tf_path, projection_path))

# Figure 6: focused AOPGraphExplorer network plus global hubs -----------------
default_network_json <- file.path(OUT_DIR, "mxenes_aop_network", "AOP_Network.json")
network_json <- path_setting("AOP_NETWORK_JSON", default_network_json)
if (file.exists(network_json)) {
  network_text <- paste(readLines(network_json, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  network_text <- gsub("(?<=:)\\s*NaN", " null", network_text, perl = TRUE)
  network_data <- jsonlite::fromJSON(network_text, simplifyDataFrame = TRUE)
  nodes <- as.data.frame(network_data$nodes, stringsAsFactors = FALSE) |>
    mutate(id = as.character(id), label = stringr::str_squish(label))
  edges <- as.data.frame(network_data$links, stringsAsFactors = FALSE) |>
    mutate(source = as.character(source), target = as.character(target))
  edges <- edges |> filter(source %in% nodes$id, target %in% nodes$id)
  vertices <- nodes |>
    transmute(
      name = id,
      label = label,
      event_type = as.character(event_type),
      aop_count = suppressWarnings(as.numeric(aop_count))
    )
  graph <- igraph::graph_from_data_frame(edges[, c("source", "target")], directed = TRUE, vertices = vertices)
  node_degree <- igraph::degree(graph, mode = "all")
  focal_ids <- intersect(c("177", "55", "188", "1392", "1496", "352", "281"), names(node_degree))
  if (!length(focal_ids)) focal_ids <- names(sort(node_degree, decreasing = TRUE))[seq_len(min(7, length(node_degree)))]
  candidates <- unique(unlist(igraph::ego(graph, order = 1, nodes = focal_ids, mode = "all")))
  candidate_ids <- igraph::V(graph)[candidates]$name
  candidate_ids <- unique(c(focal_ids, names(sort(node_degree[candidate_ids], decreasing = TRUE))[seq_len(min(42, length(candidate_ids)))]))
  focused <- igraph::induced_subgraph(graph, vids = candidate_ids)
  igraph::V(focused)$degree <- node_degree[igraph::V(focused)$name]
  igraph::V(focused)$is_focal <- igraph::V(focused)$name %in% focal_ids
  igraph::V(focused)$plot_label <- ifelse(
    igraph::V(focused)$is_focal | igraph::V(focused)$degree >= stats::quantile(igraph::V(focused)$degree, 0.85),
    stringr::str_wrap(igraph::V(focused)$label, 24), ""
  )
  set.seed(20260825)
  network_panel <- ggraph(focused, layout = "fr") +
    geom_edge_link(
      colour = "#95A1A8", alpha = 0.45, linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed"),
      end_cap = circle(2.5, "mm")
    ) +
    geom_node_point(aes(size = degree, fill = event_type, shape = is_focal), colour = "#25323A", stroke = 0.35) +
    geom_node_text(aes(label = plot_label), repel = TRUE, size = 3.1, lineheight = 0.9, max.overlaps = Inf) +
    scale_fill_manual(values = c(
      MolecularInitiatingEvent = "#7B61A8", KeyEvent = palette[["teal"]], AdverseOutcome = palette[["orange"]]
    ), labels = c(
      MolecularInitiatingEvent = "MIE", KeyEvent = "KE", AdverseOutcome = "AO"
    )) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 21), guide = "none") +
    scale_size_continuous(range = c(2.5, 8), guide = "none") +
    labs(title = "A  Focused one-hop subnetwork around recurrent hubs", fill = "Event type", size = "Degree") +
    theme_graph(base_family = "Arial") +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")

  hubs <- nodes |>
    mutate(degree = as.numeric(node_degree[id])) |>
    arrange(desc(degree)) |>
    slice_head(n = 10) |>
    mutate(label_short = stringr::str_wrap(label, 27), label_short = factor(label_short, levels = rev(label_short)))
  hub_panel <- ggplot(hubs, aes(x = degree, y = label_short, fill = event_type)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = degree), hjust = -0.22, size = 3.1, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    scale_fill_manual(values = c(
      MolecularInitiatingEvent = "#7B61A8", KeyEvent = palette[["teal"]], AdverseOutcome = palette[["orange"]]
    )) +
    labs(x = "Total degree", y = NULL, title = "B  Highest-degree events in the full export") +
    theme_publication(9.5) +
    theme(panel.grid.major.y = element_blank(), legend.position = "none", axis.text.y = element_text(size = 8.4))
  p6 <- network_panel + hub_panel + plot_layout(widths = c(1.35, 0.9)) +
    plot_annotation(
      title = sprintf("AOPGraphExplorer context: %d events and %d directed relationships", nrow(nodes), nrow(edges)),
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
    )
  save_publication_plot(p6, "Figure6_AOPGraphExplorer_readable", 13, 7.2, network_json)
}

# Figure 7: explicitly exploratory WGCNA/AOPxLINK panels ----------------------
run_root <- path_setting("AOP_AOPXGENENET_RUN_ROOT", file.path(OUT_DIR, "AOPxGeneNet_run_data"))
run_dirs <- if (dir.exists(run_root)) list.dirs(run_root, full.names = TRUE, recursive = FALSE) else character()
run_dir <- trimws(Sys.getenv("AOP_AOPXGENENET_RUN_DIR", unset = ""))
if (!nzchar(run_dir) && length(run_dirs)) run_dir <- sort(run_dirs, decreasing = TRUE)[[1]]
if (nzchar(run_dir) && dir.exists(run_dir)) {
  trait_path <- file.path(run_dir, "module_trait_cor.csv")
  trait_p_path <- file.path(run_dir, "module_trait_pvalue.csv")
  module_aop_path <- file.path(run_dir, "AOPxLink_module_AOP_enrichment.csv")
  genes_path <- file.path(run_dir, "AOPxLink_step3_selected_genes.csv")
  step2_all_stats_path <- file.path(run_dir, "AOPxLink_step2_stats_all_contrasts.csv")
  step3_all_summary_path <- file.path(run_dir, "AOPxLink_step3_all_contrasts_summary.csv")
  contrast_manifest_path <- file.path(run_dir, "contrasts", "contrasts_manifest.json")
  report_path <- file.path(run_dir, "analysis_report.json")
  required <- c(
    trait_path, trait_p_path, module_aop_path, genes_path,
    step2_all_stats_path, step3_all_summary_path, contrast_manifest_path, report_path
  )
  if (all(file.exists(required))) {
    contrast_manifest <- jsonlite::fromJSON(contrast_manifest_path, simplifyVector = FALSE)
    contrast_names <- vapply(contrast_manifest$contrasts, `[[`, character(1), "contrast")
    contrast_names <- contrast_names[contrast_names %in% paste0(CHEMS, "_vs_CTRL")]
    contrast_deg_paths <- file.path(run_dir, "deg", paste0("condition__", contrast_names), "deg_results.csv")
    if (!all(file.exists(contrast_deg_paths))) {
      stop("AOPxLINK contrast manifest refers to missing Step-1 DEG outputs.", call. = FALSE)
    }

    contrast_deg <- lapply(seq_along(contrast_names), function(i) {
      contrast_name <- contrast_names[[i]]
      must_read_csv(contrast_deg_paths[[i]]) |>
        transmute(
          Contrast = contrast_name,
          Condition = sub("_vs_CTRL$", "", contrast_name),
          gene_symbol = as.character(gene_symbol),
          log2FC = dplyr::coalesce(as.numeric(log2FoldChange), as.numeric(logFC)),
          padj = as.numeric(padj)
        )
    })
    contrast_deg <- bind_rows(contrast_deg) |>
      mutate(Condition = factor(Condition, levels = intersect(CHEMS, unique(Condition))))

    trait <- read.csv(trait_path, check.names = FALSE, row.names = 1)
    trait_p <- read.csv(trait_p_path, check.names = FALSE, row.names = 1)
    condition_cols <- setdiff(grep("^condition_", names(trait), value = TRUE), grep("^condition_factor_", names(trait), value = TRUE))
    condition_cols <- setdiff(condition_cols, "condition_CTRL")
    trait_long <- as.data.frame(trait[, condition_cols, drop = FALSE]) |>
      tibble::rownames_to_column("Module") |>
      pivot_longer(-Module, names_to = "Condition", values_to = "Correlation")
    p_long <- as.data.frame(trait_p[, condition_cols, drop = FALSE]) |>
      tibble::rownames_to_column("Module") |>
      pivot_longer(-Module, names_to = "Condition", values_to = "nominal_p")
    trait_long <- trait_long |>
      left_join(p_long, by = c("Module", "Condition")) |>
      mutate(
        Condition = sub("^condition_", "", Condition),
        Condition = factor(Condition, levels = intersect(CHEMS, unique(Condition))),
        Module = factor(Module, levels = rev(rownames(trait))),
        label = sprintf("%.2f%s", Correlation, ifelse(nominal_p < 0.001, "***", ifelse(nominal_p < 0.01, "**", ifelse(nominal_p < 0.05, "*", ""))))
      )
    trait_panel <- ggplot(trait_long, aes(x = Condition, y = Module, fill = Correlation)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = label), size = 3.0, fontface = "bold") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, fill = "Pearson r", title = "A  Module-trait correlations (nominal P)") +
      theme_publication(9.5) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

    module_aop <- must_read_csv(module_aop_path) |>
      filter(padj < FDR_CUTOFF) |>
      mutate(pair = sprintf("%s - AOP%s", module_id, aop_id), pair = factor(pair, levels = rev(pair)))
    aop_panel <- ggplot(module_aop, aes(x = log2_enrichment_ratio, y = pair)) +
      geom_segment(aes(x = 0, xend = log2_enrichment_ratio, yend = pair), colour = "#9AA3A8", linewidth = 0.8) +
      geom_point(aes(size = overlap_genes_k, fill = minus_log10_fdr), shape = 21, colour = "#25323A") +
      scale_fill_viridis_c(option = "C") +
      scale_size_continuous(range = c(3, 8)) +
      labs(x = "log2 enrichment ratio", y = NULL, size = "Overlap", fill = "-log10(FDR)", title = "B  FDR-significant module-AOP pairs") +
      theme_publication(9.5) + theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

    priority_genes <- must_read_csv(genes_path) |>
      arrange(desc(visual_priority)) |>
      distinct(gene_symbol, .keep_all = TRUE)
    gene_effects <- contrast_deg |>
      inner_join(
        priority_genes |> select(gene_symbol, module_id, visual_priority, aopxlink_score_v2),
        by = "gene_symbol"
      ) |>
      mutate(significant = padj < FDR_CUTOFF & abs(log2FC) >= LFC_MIN)
    effects_path <- file.path(PUB_DIR, "Figure7_AOPxLINK_priority_gene_effects.csv")
    readr::write_csv(gene_effects, effects_path)

    step2_all_stats <- must_read_csv(step2_all_stats_path)
    step3_all_summary <- must_read_csv(step3_all_summary_path)
    metric_value <- function(metric_name) {
      values <- step2_all_stats |>
        filter(metric == metric_name) |>
        pull(value) |>
        as.numeric()
      if (!length(values)) return(NA_real_)
      values[[1]]
    }
    semgraph_summary_paths <- file.path(
      run_dir, "contrasts", contrast_names, "step3", "semgraph", "semgraph_fit_summary.json"
    )
    semgraph_summary <- lapply(seq_along(contrast_names), function(i) {
      sem_path <- semgraph_summary_paths[[i]]
      if (!file.exists(sem_path)) {
        return(data.frame(
          Contrast = contrast_names[[i]], semgraph_fit = FALSE,
          semgraph_samples = NA_integer_, semgraph_nodes = NA_integer_, semgraph_edges = NA_integer_
        ))
      }
      sem <- jsonlite::fromJSON(sem_path)
      data.frame(
        Contrast = contrast_names[[i]], semgraph_fit = TRUE,
        semgraph_samples = sem$n_input_samples,
        semgraph_nodes = sem$n_graph_nodes,
        semgraph_edges = sem$n_graph_edges
      )
    }) |>
      bind_rows()
    whole_run_summary <- contrast_deg |>
      group_by(Contrast, Condition) |>
      summarise(
        step1_significant_genes = sum(padj < FDR_CUTOFF & abs(log2FC) >= LFC_MIN, na.rm = TRUE),
        .groups = "drop"
      ) |>
      left_join(
        gene_effects |>
          group_by(Contrast) |>
          summarise(
            prioritized_genes_significant = sum(significant, na.rm = TRUE),
            .groups = "drop"
          ),
        by = "Contrast"
      ) |>
      left_join(
        step3_all_summary |>
          transmute(
            Contrast = contrast,
            selected_genes = as.integer(selected_genes),
            selected_overlap_rows = as.integer(selected_overlap_rows),
            integrated_graph_nodes = as.integer(graph_nodes),
            integrated_graph_edges = as.integer(graph_edges)
          ),
        by = "Contrast"
      ) |>
      left_join(semgraph_summary, by = "Contrast") |>
      mutate(
        shared_wgcna_genes = as.integer(metric_value("n_gcn_nodes")),
        shared_mapped_genes = as.integer(metric_value("n_mapped_genes")),
        shared_key_events = as.integer(metric_value("n_kes")),
        shared_aops = as.integer(metric_value("n_aops")),
        shared_module_aop_pairs_tested = nrow(must_read_csv(module_aop_path)),
        shared_module_aop_pairs_fdr_lt_0_05 = sum(
          must_read_csv(module_aop_path)$padj < FDR_CUTOFF, na.rm = TRUE
        )
      )
    whole_run_summary_path <- file.path(PUB_DIR, "AOPxLINK_whole_run_summary.csv")
    readr::write_csv(whole_run_summary, whole_run_summary_path)

    top_gene_names <- priority_genes |> slice_head(n = 12) |> pull(gene_symbol)
    top_gene_effects <- gene_effects |>
      filter(gene_symbol %in% top_gene_names) |>
      mutate(gene_symbol = factor(gene_symbol, levels = rev(top_gene_names)))
    effect_limit <- max(1, max(abs(top_gene_effects$log2FC), na.rm = TRUE))
    gene_panel <- ggplot(top_gene_effects, aes(x = Condition, y = gene_symbol, fill = log2FC)) +
      geom_tile(colour = "white", linewidth = 0.55) +
      geom_point(
        data = top_gene_effects |> filter(significant),
        shape = 8, size = 2.2, colour = "#17232B"
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-effect_limit, effect_limit), oob = scales::squish
      ) +
      labs(
        x = NULL, y = NULL, fill = "Step-1\nlog2 FC",
        title = "C  Run-level prioritized genes across contrasts",
        subtitle = sprintf("Asterisk: limma FDR < %.2f and |log2 FC| >= %.1f", FDR_CUTOFF, LFC_MIN)
      ) +
      theme_publication(9.2) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom"
      )

    p7 <- trait_panel + aop_panel + gene_panel +
      plot_layout(widths = c(1.15, 0.85, 1.1), guides = "keep") +
      plot_annotation(
        title = "Integrated five-contrast AOPxGeneNet/AOPxLINK analysis",
        subtitle = "Whole-run synthesis of shared coexpression/AOP structure and contrast-specific Step-1 effects",
        theme = theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 10))
      ) & theme(legend.position = "bottom")
    figure7_sources <- c(required, contrast_deg_paths, semgraph_summary_paths, effects_path, whole_run_summary_path)
    save_publication_plot(p7, "Figure7_AOPxGeneNet_exploratory_readable", 13.5, 7.5, figure7_sources)

    # Supplementary figure: display the complete 31-gene run-level priority set
    # against the condition-specific Step-1 effects for all five contrasts.
    full_gene_order <- priority_genes$gene_symbol
    full_gene_effects <- gene_effects |>
      mutate(
        gene_label = sprintf("%s (%s)", gene_symbol, module_id),
        gene_label = factor(
          gene_label,
          levels = rev(sprintf(
            "%s (%s)", full_gene_order,
            priority_genes$module_id[match(full_gene_order, priority_genes$gene_symbol)]
          ))
        )
      )
    full_effect_limit <- max(1, max(abs(full_gene_effects$log2FC), na.rm = TRUE))
    p_s2 <- ggplot(full_gene_effects, aes(x = Condition, y = gene_label, fill = log2FC)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_point(
        data = full_gene_effects |> filter(significant),
        shape = 8, size = 2.3, colour = "#17232B"
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-full_effect_limit, full_effect_limit), oob = scales::squish
      ) +
      labs(
        x = NULL, y = NULL, fill = "Step-1 log2 FC",
        title = "All run-level AOPxLINK-prioritized genes across five contrasts",
        subtitle = sprintf(
          "Rows show gene (WGCNA module); asterisks mark limma FDR < %.2f and |log2 FC| >= %.1f",
          FDR_CUTOFF, LFC_MIN
        )
      ) +
      theme_publication(10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8.3),
        legend.position = "bottom"
      )
    save_publication_plot(
      p_s2, "FigureS2_AOPxLINK_prioritized_gene_effects", 8.5, 10.5,
      c(genes_path, contrast_deg_paths, effects_path, whole_run_summary_path)
    )
  }
}

# Supplementary specificity figure generated by step 11 -----------------------
specificity_dir <- file.path(OUT_DIR, "Specificity_null")
specificity_png <- file.path(specificity_dir, "Figure_AOP_specificity_null.png")
specificity_pdf <- file.path(specificity_dir, "Figure_AOP_specificity_null.pdf")
if (file.exists(specificity_png) && file.exists(specificity_pdf)) {
  target_png <- file.path(PUB_DIR, "FigureS1_AOP_specificity_null.png")
  target_pdf <- file.path(PUB_DIR, "FigureS1_AOP_specificity_null.pdf")
  file.copy(specificity_png, target_png, overwrite = TRUE)
  file.copy(specificity_pdf, target_pdf, overwrite = TRUE)
  manifest[[length(manifest) + 1L]] <- data.frame(
    figure = "FigureS1_AOP_specificity_null",
    png = normalizePath(target_png, winslash = "/", mustWork = TRUE),
    pdf = normalizePath(target_pdf, winslash = "/", mustWork = TRUE),
    source_files = normalizePath(file.path(specificity_dir, "specificity_summary.csv"), winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest_df <- dplyr::bind_rows(manifest)
if (!nrow(manifest_df)) stop("No publication figures were generated.", call. = FALSE)
readr::write_csv(manifest_df, file.path(PUB_DIR, "publication_figure_manifest.csv"))

cat("Publication figures complete:", normalizePath(PUB_DIR, winslash = "/"), "\n")

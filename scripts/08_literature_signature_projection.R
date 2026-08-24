#!/usr/bin/env Rscript
# Project literature-derived UP/DOWN signatures onto DESeq2 contrasts with
# preranked FGSEA, then export normalized enrichment scores and a heatmap.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "AnnotationDbi", "org.Hs.eg.db", "fgsea", "readr", "dplyr", "tidyr",
  "tibble", "pheatmap"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

PROJECTION_DIR <- file.path(OUT_DIR, "Signature_Projection")
dir.create(PROJECTION_DIR, recursive = TRUE, showWarnings = FALSE)

de_files <- setNames(
  file.path(DE_DIR, sprintf("DEG_%s_vs_CTRL_all.csv", CHEMS)),
  CHEMS
)
assert_files(c(LIT_FILE, de_files), "Signature-projection input files")

make_named_map <- function(data, key_col, value_col) {
  keep <- !is.na(data[[key_col]]) & nzchar(data[[key_col]]) &
    !is.na(data[[value_col]]) & nzchar(data[[value_col]])
  data <- data[keep, c(key_col, value_col), drop = FALSE]
  data <- data[order(data[[key_col]], data[[value_col]]), , drop = FALSE]
  data <- data[!duplicated(data[[key_col]]), , drop = FALSE]
  setNames(data[[value_col]], data[[key_col]])
}

quiet_stdout <- function(expr) {
  result <- NULL
  invisible(utils::capture.output(result <- expr, type = "output"))
  result
}

normalize_direction <- function(x) {
  text <- tolower(trimws(as.character(x)))
  numeric_value <- suppressWarnings(as.numeric(text))
  dplyr::case_when(
    is.finite(numeric_value) & numeric_value > 0 ~ "UP",
    is.finite(numeric_value) & numeric_value < 0 ~ "DOWN",
    text %in% c("up", "upregulated", "up-regulated", "increase", "increased") ~ "UP",
    text %in% c("down", "downregulated", "down-regulated", "decrease", "decreased") ~ "DOWN",
    TRUE ~ NA_character_
  )
}

lit <- suppressMessages(readr::read_csv(LIT_FILE, show_col_types = FALSE))
required_lit_columns <- c("Chemical", "Direction")
missing_lit_columns <- setdiff(required_lit_columns, names(lit))
if (length(missing_lit_columns)) {
  stop(
    "Literature table is missing columns required by step 08: ",
    paste(missing_lit_columns, collapse = ", "),
    call. = FALSE
  )
}

if ("Organism" %in% names(lit)) {
  lit <- lit %>%
    filter(grepl("human|homo sapiens|h\\. sapiens", Organism, ignore.case = TRUE))
}
lit$Direction_normalized <- normalize_direction(lit$Direction)
lit <- lit %>% filter(!is.na(Direction_normalized))

gene_col <- intersect(
  c("Gene", "SYMBOL", "symbol", "GeneSymbol", "gene", "Feature"),
  names(lit)
)
gene_col <- if (length(gene_col)) gene_col[[1]] else NULL
if (is.null(gene_col) && !"ENSG" %in% names(lit)) {
  stop("Literature table needs a gene-symbol column or `ENSG`.", call. = FALSE)
}

ensg_input <- if ("ENSG" %in% names(lit)) {
  strip_v(trimws(as.character(lit$ENSG)))
} else {
  rep(NA_character_, nrow(lit))
}
symbol_input <- if (!is.null(gene_col)) {
  toupper(trimws(as.character(lit[[gene_col]])))
} else {
  rep(NA_character_, nrow(lit))
}

valid_ensembl <- intersect(
  unique(ensg_input[!is.na(ensg_input) & nzchar(ensg_input)]),
  AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
)
ens_map_table <- if (length(valid_ensembl)) {
  suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_ensembl,
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  ))
} else {
  data.frame(ENSEMBL = character(), SYMBOL = character())
}
ens_to_symbol <- make_named_map(ens_map_table, "ENSEMBL", "SYMBOL")

official_symbols <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
valid_aliases <- intersect(
  unique(symbol_input[!is.na(symbol_input) & nzchar(symbol_input)]),
  AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ALIAS")
)
alias_map_table <- if (length(valid_aliases)) {
  suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_aliases,
    keytype = "ALIAS",
    columns = "SYMBOL"
  ))
} else {
  data.frame(ALIAS = character(), SYMBOL = character())
}
alias_to_symbol <- make_named_map(alias_map_table, "ALIAS", "SYMBOL")

symbol_from_ensembl <- unname(ens_to_symbol[ensg_input])
symbol_direct <- ifelse(symbol_input %in% official_symbols, symbol_input, NA_character_)
symbol_from_alias <- unname(alias_to_symbol[symbol_input])
lit$SYMBOL <- dplyr::coalesce(
  symbol_from_ensembl,
  symbol_direct,
  symbol_from_alias
)

literature_membership <- lit %>%
  filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
  transmute(
    Signature = paste(Chemical, Direction_normalized, sep = "_"),
    Chemical = as.character(Chemical),
    Direction = Direction_normalized,
    SYMBOL = toupper(SYMBOL)
  ) %>%
  distinct()

if (!nrow(literature_membership)) {
  stop("No human literature UP/DOWN signatures could be constructed.", call. = FALSE)
}

readr::write_csv(
  literature_membership,
  file.path(PROJECTION_DIR, "literature_signature_membership.csv")
)
literature_sets <- split(literature_membership$SYMBOL, literature_membership$Signature)
literature_sets <- lapply(literature_sets, unique)

de_tables <- lapply(de_files, function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
})

required_de_columns <- c("ENSEMBL", "log2FoldChange", "pvalue")
for (chemical in names(de_tables)) {
  missing_de_columns <- setdiff(required_de_columns, names(de_tables[[chemical]]))
  if (length(missing_de_columns)) {
    stop(
      basename(de_files[[chemical]]), " is missing: ",
      paste(missing_de_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

de_ensembl <- unique(unlist(lapply(de_tables, function(data) strip_v(data$ENSEMBL))))
de_ensembl <- de_ensembl[!is.na(de_ensembl) & nzchar(de_ensembl)]
valid_de_ensembl <- intersect(
  de_ensembl,
  AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
)
de_map_table <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = valid_de_ensembl,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
))
de_to_symbol <- make_named_map(de_map_table, "ENSEMBL", "SYMBOL")

make_rank <- function(data) {
  ensembl <- strip_v(data$ENSEMBL)
  symbol <- unname(de_to_symbol[ensembl])
  fallback_stat <- sign(data$log2FoldChange) *
    -log10(pmax(data$pvalue, .Machine$double.xmin))
  rank_stat <- if ("stat" %in% names(data)) {
    ifelse(is.finite(data$stat), data$stat, fallback_stat)
  } else {
    fallback_stat
  }

  ranked <- tibble::tibble(SYMBOL = symbol, rank_stat = rank_stat) %>%
    filter(!is.na(SYMBOL), nzchar(SYMBOL), is.finite(rank_stat)) %>%
    group_by(SYMBOL) %>%
    slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rank_stat))

  stats <- ranked$rank_stat
  names(stats) <- ranked$SYMBOL
  stats
}

ranks <- lapply(de_tables, make_rank)
message(
  "Rank sizes: ",
  paste(sprintf("%s=%d", names(ranks), lengths(ranks)), collapse = ", ")
)

flatten_list_columns <- function(data) {
  data <- as.data.frame(data)
  for (column in names(data)) {
    if (is.list(data[[column]])) {
      data[[column]] <- vapply(
        data[[column]],
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      )
    }
  }
  data
}

fgsea_results <- lapply(names(ranks), function(chemical) {
  result <- quiet_stdout(
    fgsea::fgseaMultilevel(
      pathways = literature_sets,
      stats = ranks[[chemical]],
      minSize = 10,
      maxSize = 1000,
      eps = 0,
      nproc = 1
    )
  )
  result$Contrast <- chemical
  result
})
names(fgsea_results) <- names(ranks)

for (chemical in names(fgsea_results)) {
  readr::write_csv(
    flatten_list_columns(fgsea_results[[chemical]]),
    file.path(PROJECTION_DIR, sprintf("fgsea_%s_vs_CTRL_literature_sets.csv", chemical))
  )
}

combined <- dplyr::bind_rows(fgsea_results)
if (!nrow(combined)) {
  stop("FGSEA returned no testable literature signatures.", call. = FALSE)
}
readr::write_csv(
  flatten_list_columns(combined),
  file.path(PROJECTION_DIR, "fgsea_all_contrasts_literature_sets.csv")
)

nes_wide <- combined %>%
  select(pathway, Contrast, NES) %>%
  distinct() %>%
  pivot_wider(names_from = Contrast, values_from = NES)
readr::write_csv(
  nes_wide,
  file.path(PROJECTION_DIR, "Literature_vs_GEO_NES_matrix.csv")
)

nes_matrix <- nes_wide %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()
heatmap_colors <- grDevices::colorRampPalette(c("steelblue", "white", "firebrick"))(101)
heatmap_breaks <- seq(-max(abs(nes_matrix), na.rm = TRUE),
                      max(abs(nes_matrix), na.rm = TRUE),
                      length.out = 102)

pheatmap::pheatmap(
  nes_matrix,
  main = "Literature signatures projected onto transcriptomic contrasts (NES)",
  color = heatmap_colors,
  breaks = heatmap_breaks,
  filename = file.path(PROJECTION_DIR, "Literature_vs_GEO_NES_heatmap.png"),
  width = 8,
  height = 6
)
pheatmap::pheatmap(
  nes_matrix,
  main = "Literature signatures projected onto transcriptomic contrasts (NES)",
  color = heatmap_colors,
  breaks = heatmap_breaks,
  filename = file.path(PROJECTION_DIR, "Literature_vs_GEO_NES_heatmap.pdf"),
  width = 8,
  height = 6
)

message("Signature projection complete: ", normalizePath(PROJECTION_DIR, winslash = "/"))

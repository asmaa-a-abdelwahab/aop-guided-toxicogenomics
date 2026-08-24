#!/usr/bin/env Rscript
# =============================================================================
# TF expansion from a CSV of genes (Ensembl or SYMBOL) using DoRothEA regulons
# Outputs are written under outputs/literature_tf_expansion/<run_name>/
# Default input: data/Literature_DEGs.csv
# =============================================================================

suppressPackageStartupMessages({
  source(file.path("scripts", "_config.R"))
  assert_packages(c(
    "optparse", "readr", "dplyr", "tibble", "stringr", "ggplot2",
    "AnnotationDbi", "org.Hs.eg.db", "jsonlite", "fs", "dorothea",
    "forcats", "openxlsx"
  ))
  library(optparse)
})

# ---------- CLI ----------------------------------------------------------------
opt_list <- list(
  make_option(c("-i", "--input"),      type="character", default=LIT_FILE,
              help="Input CSV path [default: %default]"),
  make_option(c("-c", "--gene-col"),   type="character", default=NA,
              help="Name of gene column to use (auto-detect if omitted)"),
  make_option(c("-t", "--id-type"),    type="character", default=NA,
              help="Input ID type: 'ensembl' or 'symbol' (auto-detect if omitted)"),
  make_option(c("-o", "--outdir"),     type="character", default=file.path(OUT_DIR, "literature_tf_expansion"),
              help="Output directory root [default: %default]"),
  make_option(c("-n", "--name"),       type="character", default=NA,
              help="Run name for the output folder (defaults to input file stem)"),
  make_option(c("--no-plot"),          action="store_false", dest="plot", default=TRUE,
              help="Skip summary bar plots"),
  make_option(c("--quiet"),            action="store_true", default=FALSE,
              help="Reduce console output")
)
opt <- parse_args(OptionParser(option_list = opt_list))

msg <- function(...) if (!isTRUE(opt$quiet)) message(...)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
  library(ggplot2); library(AnnotationDbi); library(org.Hs.eg.db)
  library(jsonlite); library(fs); library(dorothea); library(forcats)
})


# ---------- Helpers -------------------------------------------------------------
auto_detect_gene_col <- function(df){
  cands <- c("Gene","Feature","SYMBOL","symbol","gene","GeneSymbol",
             "Ensembl","ENSEMBL","EnsemblID","ENSEMBL_ID","ID",
             "gene_id","GeneID","Ensembl.Gene.ID","Ensembl_ID","ids")
  ok <- intersect(cands, names(df))
  if (length(ok)) return(ok[[1]])
  # fallback: first column
  names(df)[1]
}

clean_ens <- function(x){
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)                 # drop version suffix .xx
  x[grepl("^ENSG\\d+", x, ignore.case = FALSE)]
}

detect_id_type <- function(x){
  x <- na.omit(trimws(as.character(x)))
  frac_ens <- mean(grepl("^ENSG\\d+", x))
  if (is.finite(frac_ens) && frac_ens >= 0.7) "ensembl" else "symbol"
}

map_ensembl_to_symbol <- function(ensg){
  if (!length(ensg)) return(setNames(character(0), character(0)))
  ensg <- unique(ensg)
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=ensg, keytype="ENSEMBL", columns=c("SYMBOL"))
  m <- m[!is.na(m$SYMBOL) & nzchar(m$SYMBOL), ]
  m <- m[!duplicated(m$ENSEMBL), c("ENSEMBL","SYMBOL")]
  setNames(m$SYMBOL, m$ENSEMBL)
}

map_symbol_to_ensembl <- function(sym){
  if (!length(sym)) return(setNames(character(0), character(0)))
  sym <- unique(sym)
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=sym, keytype="SYMBOL", columns=c("ENSEMBL"))
  m <- m[!is.na(m$ENSEMBL) & nzchar(m$ENSEMBL), ]
  m <- m[!duplicated(m$SYMBOL), c("SYMBOL","ENSEMBL")]
  setNames(m$ENSEMBL, m$SYMBOL)
}

# Expand input genes (Ensembl or SYMBOL) with TF targets (DoRothEA)
expand_with_tf_targets <- function(genes_in,
                                   id_type = c("ensembl","symbol"),
                                   regulons,
                                   return_id_type = c("ensembl","symbol")){
  id_type <- match.arg(id_type)
  return_id_type <- match.arg(return_id_type)

  genes_in <- unique(na.omit(trimws(as.character(genes_in))))
  if (!length(genes_in)) {
    return(list(
      expanded      = tibble(Gene = character()),
      original      = tibble(Gene = character()),
      tf_hits       = regulons[0,],
      targets_only  = tibble(Gene = character()),
      tf_summary    = tibble(),
      note          = "Empty input"
    ))
  }

  # Work in SYMBOL-space for matching regulons
  if (id_type == "ensembl"){
    ens2sym <- map_ensembl_to_symbol(genes_in)
    symbols_in <- unique(na.omit(unname(ens2sym[genes_in])))
  } else {
    symbols_in <- unique(toupper(genes_in))
  }

  tf_hits <- regulons %>% dplyr::filter(toupper(tf) %in% symbols_in)
  if (nrow(tf_hits) == 0){
    # No TFs -> no expansion, just echo input
    out <- if (return_id_type == "ensembl") genes_in else unique(na.omit(unname(map_ensembl_to_symbol(genes_in))))
    return(list(
      expanded      = tibble(Gene = sort(out)),
      original      = tibble(Gene = sort(if (return_id_type=="ensembl") genes_in else toupper(symbols_in))),
      tf_hits       = regulons[0,],
      targets_only  = tibble(Gene = character()),
      tf_summary    = tibble(),
      note          = "No TFs matched DoRothEA at chosen confidence."
    ))
  }

  targets_sym <- sort(unique(toupper(tf_hits$target)))
  added_sym   <- sort(setdiff(targets_sym, symbols_in))
  all_sym     <- sort(unique(c(symbols_in, targets_sym)))

  if (return_id_type == "ensembl"){
    sym2ens <- map_symbol_to_ensembl(all_sym)
    expanded <- unique(na.omit(unname(sym2ens[all_sym])))
    sym2ens_targets <- map_symbol_to_ensembl(added_sym)
    targets_only <- unique(na.omit(unname(sym2ens_targets[added_sym])))
    original <- if (id_type == "ensembl") genes_in else unique(na.omit(unname(map_symbol_to_ensembl(symbols_in))))
  } else {
    expanded <- all_sym
    targets_only <- added_sym
    original <- if (id_type == "symbol") symbols_in else unique(na.omit(unname(map_ensembl_to_symbol(genes_in))))
  }

  tf_summary <- tf_hits %>%
    mutate(tf = toupper(tf), target = toupper(target)) %>%
    group_by(tf) %>%
    summarise(n_targets = n_distinct(target),
              n_added   = n_distinct(setdiff(target, symbols_in)),
              confidences = paste(sort(unique(confidence)), collapse = ","),
              targets_added = paste(sort(setdiff(target, symbols_in)), collapse = ";"),
              .groups = "drop") %>%
    arrange(desc(n_added), desc(n_targets))

  list(
    expanded      = tibble(Gene = expanded),
    original      = tibble(Gene = sort(unique(original))),
    tf_hits       = tf_hits %>% arrange(tf, target),
    targets_only  = tibble(Gene = targets_only),
    tf_summary    = tf_summary,
    note          = paste0("Expanded with DoRothEA in SYMBOL-space; returned as ", toupper(return_id_type), ".")
  )
}

# ---------- Load input ----------------------------------------------------------
assert_files(opt$input, "Literature DEG input")
msg("Reading: ", opt$input)
df_in <- suppressMessages(readr::read_csv(opt$input, show_col_types = FALSE))

chem_col <- "Chemical"  # or detect it
if (!chem_col %in% names(df_in)) {
  stop("Input must contain a 'Chemical' column.", call. = FALSE)
}
chemicals <- df_in[[chem_col]] |> as.character() |> trimws()
chemicals <- chemicals[nzchar(chemicals)]              # drop empty strings
chemicals <- unique(chemicals)
chemicals <- sort(chemicals)
if (!length(chemicals)) stop("No non-empty chemicals found in the input.", call. = FALSE)


for (chem in chemicals) {
  chem_col <- "Chemical"  # or detect it if needed
  df_chem <- df_in[ !is.na(df_in[[chem_col]]) &
                      toupper(trimws(as.character(df_in[[chem_col]]))) ==
                      toupper(trimws(as.character(chem))), ,
                    drop = FALSE ]

  gene_col <- if (is.na(opt$`gene-col`)) auto_detect_gene_col(df_chem) else opt$`gene-col`
  if (!gene_col %in% names(df_chem)) stop("Gene column not found: ", gene_col)

  genes_raw <- df_chem[[gene_col]]

  # Guess ID type if not specified
  id_type <- if (is.na(opt$`id-type`)) detect_id_type(genes_raw) else tolower(opt$`id-type`)
  if (!id_type %in% c("ensembl","symbol")) stop("id-type must be 'ensembl' or 'symbol'")

  if (id_type == "ensembl") {
    genes <- clean_ens(genes_raw)
  } else {
    genes <- toupper(trimws(as.character(genes_raw)))
    genes <- genes[nzchar(genes) & !is.na(genes)]
  }

  # ---------- Output location per run -----------------------------------------
  run_name <- if (!is.na(opt$name)) opt$name else tools::file_path_sans_ext(basename(opt$input))
  out_dir  <- fs::path(opt$outdir, run_name)
  fs::dir_create(out_dir, recurse = TRUE)
  
  # ---------- Run both confidence bundles -------------------------------------
  bundles <- list(AB = c("A","B"),
                  ABC = c("A","B","C"))
  
  for (tag in names(bundles)) {
    keep_conf <- bundles[[tag]]
    
    # ---------- Regulons subset (schema-robust) --------------------------------
    reg <- dorothea::dorothea_hs
    if (!"likelihood" %in% names(reg)) reg$likelihood <- NA_real_
    if (!"mor" %in% names(reg)) {
      reg$mor <- if ("mor_sign" %in% names(reg)) reg$mor_sign else NA_character_
    }
    if ("confidence" %in% names(reg)) reg$confidence <- toupper(reg$confidence)
    
    dorothea_regulons <- reg %>%
      dplyr::filter(confidence %in% toupper(keep_conf)) %>%
      dplyr::select(dplyr::any_of(c("tf","target","mor","likelihood","confidence")))
    
    msg(sprintf("Regulons: %s edges at confidence {%s} for %s",
                nrow(dorothea_regulons), paste(keep_conf, collapse=","), tag))
    
    # ---------- Expansion -------------------------------------------------------
    exp_res <- expand_with_tf_targets(
      genes_in        = genes,
      id_type         = id_type,
      regulons        = dorothea_regulons,
      return_id_type  = "ensembl"
    )
    
    # ---------- Save outputs (suffix by bundle) --------------------------------
    readr::write_csv(exp_res$original,     fs::path(out_dir, sprintf("TFEXP_%s_ALL_ORIG_ORIGINAL.csv", chem)))
    readr::write_csv(exp_res$expanded,     fs::path(out_dir, sprintf("TFEXP_%s_ALL_%s_EXPANDED.csv", chem, tag)))
    readr::write_csv(exp_res$targets_only, fs::path(out_dir, sprintf("TFEXP_%s_ALL_%s_TARGETS_ONLY.csv", chem, tag)))
    readr::write_csv(exp_res$tf_hits,      fs::path(out_dir, sprintf("TFEXP_%s_ALL_%s_HITS.csv",    chem, tag)))
    readr::write_csv(exp_res$tf_summary,   fs::path(out_dir, sprintf("TFEXP_%s_ALL_%s_SUMMARY.csv",  chem, tag)))
    
    report <- list(
      input_file   = normalizePath(opt$input),
      chemical     = chem,
      gene_column  = gene_col,
      id_type      = id_type,
      n_input      = length(unique(genes)),
      n_expanded   = nrow(exp_res$expanded),
      n_added      = nrow(exp_res$targets_only),
      n_TFs_hit    = dplyr::n_distinct(exp_res$tf_hits$tf),
      confidence   = keep_conf,
      note         = exp_res$note,
      tag          = tag,
      timestamp    = as.character(Sys.time())
    )
    
    writeLines(c(
      "=== TF Expansion Report ===",
      paste("Input:",      report$input_file),
      paste("Chemical:",   report$chemical),
      paste("Gene column:",report$gene_column),
      paste("ID type:",    report$id_type),
      paste("Input genes:",report$n_input),
      paste("Expanded genes:", report$n_expanded),
      paste("Targets added:",  report$n_added),
      paste("TFs matched:",    report$n_TFs_hit),
      paste("Confidence:",     paste(report$confidence, collapse=",")),
      paste("Bundle tag:",     report$tag),
      paste("Note:",           report$note),
      paste("Time:",           report$timestamp)
    ), con = fs::path(out_dir, sprintf("%s_expansion_report_%s.txt",  chem, tag)))
    
    write(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE),
          file = fs::path(out_dir, sprintf("%s_expansion_report_%s.json", chem, tag)))
    
    # ---------- Optional plot ---------------------------------------------------
    if (isTRUE(opt$plot) && nrow(exp_res$tf_summary)) {
      p <- exp_res$tf_summary %>%
        mutate(tf = forcats::fct_reorder(tf, n_added)) %>%
        ggplot(aes(x = n_added, y = tf)) +
        geom_col() +
        labs(
          title = sprintf("TFs ranked by NEW targets added (%s, %s)", chem, tag),
          x = "Targets added",
          y = "Transcription Factor"
        ) +
        theme_minimal(base_size = 14)
      
      ggsave(fs::path(out_dir, sprintf("%s_tf_targets_added_bar_%s.png", chem, tag)),
             p, width = 10, height = 8, dpi = 600)
      ggsave(fs::path(out_dir, sprintf("%s_tf_targets_added_bar_%s.pdf", chem, tag)),
             p, width = 10, height = 8, device = cairo_pdf)
    }
  }
}

msg("✓ Done. Files saved under: ", normalizePath(out_dir))


# ===================== Merge literature DEG sets =============================
# - Discovers all original and expanded CSVs created above
# - Adds Experiment = source filename (without extension)
# - Renames the gene column to Feature (auto-detects common variants)
# - Adds timepoint = 24
# - Writes Excel: Literature_DEGs_Expanded.xlsx
# ============================================================================
files <- list.files(
  out_dir,
  pattern = "^TFEXP_.*_(EXPANDED|ORIGINAL)\\.csv$",
  full.names = TRUE
)
if (!length(files)) stop("No TF expansion tables were produced.", call. = FALSE)

# --- helper: find & standardize the gene column to "Feature" -----------------
pick_gene_col <- function(nms) {
  # common variants (case-insensitive)
  candidates <- c(
    "feature","gene","genes","symbol","gene_symbol","hgnc_symbol",
    "external_gene_name","gene_name","geneid","gene_id",
    "Gene","Genes","Symbol","Gene.Symbol","HGNC.symbol","SYMBOL","GENE","Feature"
  )
  idx <- match(tolower(nms), tolower(candidates))
  hit <- which(!is.na(idx))[1]
  if (length(hit) == 0) NA_character_ else nms[hit]
}

read_one <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  exp_name <- tools::file_path_sans_ext(basename(path))
  df <- dplyr::mutate(df, Experiment = exp_name, .before = 1)
  gcol <- pick_gene_col(names(df))
  if (!is.na(gcol) && gcol != "Feature") {
    df <- dplyr::rename(df, Feature = all_of(gcol))
  } else if (is.na(gcol)) {
    warning(sprintf("No gene-like column found in: %s", basename(path)))
  }
  df
}

# --- merge -------------------------------------------------------------------
merged <- dplyr::bind_rows(lapply(files, read_one))

# Add the case-study timepoint configured in scripts/_config.R.
merged <- merged %>% mutate(timepoint = TIMEPOINT, .after = Experiment)

# --- write Excel if possible; otherwise CSV fallback -------------------------
write_excel <- function(df, path_xlsx) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(df, path_xlsx, asTable = TRUE)
    message("Saved Excel with openxlsx: ", normalizePath(path_xlsx))
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(df, path_xlsx)
    message("Saved Excel with writexl: ", normalizePath(path_xlsx))
  } else {
    warning("Neither 'openxlsx' nor 'writexl' is installed. Saving CSV instead.")
    path_csv <- sub("\\.xlsx$", ".csv", path_xlsx)
    readr::write_csv(df, path_csv)
    message("Saved CSV: ", normalizePath(path_csv))
  }
}


# ---- 5) Write Excel ---------------------------------------------------------
out_xlsx <- file.path(DE_DIR, "Literature_DEGs_Expanded.xlsx")
write_excel(merged, path_xlsx = out_xlsx)
message("Saved: ", normalizePath(out_xlsx))

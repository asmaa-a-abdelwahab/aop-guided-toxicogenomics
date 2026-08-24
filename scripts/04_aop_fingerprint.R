#!/usr/bin/env Rscript
# Key Event enrichment, AOP fingerprinting, and cross-chemical comparisons.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "AOPfingerprintR", "AnnotationDbi", "org.Hs.eg.db", "pheatmap",
  "ggplot2", "readr", "readxl", "dplyr", "tidyr", "tibble",
  "htmlwidgets", "htmltools", "visNetwork"
))

suppressPackageStartupMessages({
  library(AOPfingerprintR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readxl)
  library(dplyr)
})

# ----- Single-folder output (everything goes under outputs/AOP) -----
dir.create(AOP_DIR, recursive = TRUE, showWarnings = FALSE)


# reference sets
data("human_ens_clusters", package = "AOPfingerprintR")  # KEs
data("human_ens_aop",     package = "AOPfingerprintR")  # AO chains

# --- Helper to normalize ENSG IDs (strip version) ---
.norm_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  unique(x[nzchar(x) & !is.na(x)])
}

# Some AOPfingerprintR functions print internal labels directly to stdout even
# with verbose = FALSE. Capture stdout only; messages, warnings, and errors stay
# visible to the pipeline.
quiet_stdout <- function(expr) {
  result <- NULL
  invisible(utils::capture.output(result <- expr, type = "output"))
  result
}

# --- Load project-specific backgrounds ---
# background_genes.xlsx        -> for GSE chemicals
# background_genes_GRCh38.xlsx -> for literature chemicals
background_files <- c(
  file.path(ROOT, "data", "background_genes.xlsx"),
  file.path(ROOT, "data", "background_genes_GRCh38.xlsx")
)
assert_files(background_files, "Background gene files")
background <- suppressMessages(
  read_excel(background_files[[1]], col_names = FALSE)
)[[1]]
genome_background <- suppressMessages(
  read_excel(background_files[[2]], col_names = FALSE)
)[[1]]


bg_gse <- .norm_ids(background)
bg_lit <- .norm_ids(genome_background)


message(sprintf("Background sizes (AS-IS) | GSE=%d  LIT=%d",
                length(bg_gse), length(bg_lit)))


# ---- Build GList from DE outputs --------------------------------------------
to_glist <- function(vec, chem, setlab){
  if (!length(vec)) return(NULL)
  data.frame(Feature = vec, timepoint = 24L, Experiment = paste0(chem, "_", setlab))
}

GList_gse <- list()
for (chem in CHEMS) {
  for (setlab in RUN_ON) {
    path <- file.path(DE_DIR, sprintf("DEG_%s_%s_ensg.txt", chem, setlab))
    if (file.exists(path)) {
      key <- paste0(chem, "_", setlab)
      GList_gse[[key]] <- to_glist(readr::read_lines(path), chem, setlab)
    }
  }
}
GList_gse <- Filter(Negate(is.null), GList_gse)
if (!length(GList_gse)) {
  stop("No DEG gene lists found. Run step 01 before step 04.", call. = FALSE)
}

# --- Normalize Feature IDs (strip Ensembl version suffix) --------------------
normalize_glist <- function(GList) {
  lapply(GList, function(df) {
    df$Feature <- sub("\\.\\d+$", "", as.character(df$Feature))
    unique(df)
  })
}

# Save enrichment results whether a data.frame or a list of data.frames
save_enrich <- function(obj, basepath) {
  if (is.null(obj)) return(invisible(NULL))
  dir.create(dirname(basepath), showWarnings = FALSE, recursive = TRUE)

  if (is.data.frame(obj)) {
    utils::write.csv(obj, paste0(basepath, ".csv"), row.names = FALSE)
    return(obj)
  }

  # list case: bind all non-empty data.frames
  parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, obj)
  if (length(parts)) {
    df <- dplyr::bind_rows(parts, .id = "InputSet")
    utils::write.csv(df, paste0(basepath, "_BOUND.csv"), row.names = FALSE)
    return(df)
  } else {
    saveRDS(obj, paste0(basepath, ".rds"))
    return(invisible(NULL))
  }
}

GList_gse <- normalize_glist(GList_gse)

enrich_ke <- function(GList, label, out_dir = AOP_DIR, background) {
  if (!length(GList)) return(invisible(NULL))

  ke_sig <- quiet_stdout(
    suppressPackageStartupMessages(enrich_KEs_AOPs(
      GList            = GList,
      list_gene_sets   = human_ens_clusters,
      only_significant = TRUE,   # <- only significant
      pval_th          = 0.05,
      adj.method       = "fdr",
      merge_by         = "Ke",
      background       = background,    # <- use provided background AS-IS
      numerical_properties = NULL,
      verbose          = FALSE
    ))
  )
  # Save if it's a data.frame or a list of data.frames
  if (is.null(ke_sig)) return(invisible(NULL))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.data.frame(ke_sig)) {
    utils::write.csv(ke_sig, file.path(out_dir, paste0("KE_enrichment_", label, "_SIG.csv")), row.names = FALSE)
  } else {
    parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, ke_sig)
    if (length(parts)) {
      df <- dplyr::bind_rows(parts, .id = "InputSet")
      utils::write.csv(df, file.path(out_dir, paste0("KE_enrichment_", label, "_SIG_BOUND.csv")), row.names = FALSE)
    } else {
      saveRDS(ke_sig, file.path(out_dir, paste0("KE_enrichment_", label, "_SIG.rds")))
    }
  }
  invisible(ke_sig)
}

ke_gse <- enrich_ke(GList_gse, "GSE_DEGs", background = bg_gse)



# ---- Literature branch (optional if file exists) ----------------------------
LIT_FILE <- file.path(DE_DIR, "Literature_DEGs_Expanded.xlsx")

# helper: pick gene column robustly
pick_gene_col <- function(df) {
  cands <- c("Feature","Gene","SYMBOL","symbol","gene","ENSEMBL","Ensembl","GeneSymbol","GeneID","ID")
  ok <- intersect(cands, names(df))
  if (length(ok)) ok[[1]] else stop("Couldn't find a gene column. Columns present: ", paste(names(df), collapse=", "))
}

# ensure Ensembl IDs (map from SYMBOL only when needed)
to_ensembl <- function(x_chr) {
  x_chr <- trimws(as.character(x_chr))
  x_chr <- sub("\\.\\d+$", "", x_chr)                 # drop .version if present
  ens_like <- grepl("^ENSG\\d+$", x_chr)
  if (all(ens_like)) return(unique(x_chr[ens_like]))
  # SYMBOL -> ENSEMBL
  syms <- unique(x_chr[!ens_like & nzchar(x_chr)])
  if (!length(syms)) return(character(0))
  m <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                          keytype = "SYMBOL", columns = "ENSEMBL")
  )
  unique(na.omit(m$ENSEMBL))
}

if (file.exists(LIT_FILE)) {
  # read Excel
  lit <- readxl::read_excel(LIT_FILE)
  
  # filter to human if the column exists
  if ("Organism" %in% names(lit)) {
    lit <- lit %>% dplyr::filter(grepl("human|homo sapiens", Organism, ignore.case = TRUE))
  }
  
  # decide grouping column: prefer Chemical; else use Experiment
  group_col <- if ("Chemical" %in% names(lit)) "Chemical"
  else if ("Experiment" %in% names(lit)) "Experiment"
  else stop("Expected column 'Chemical' or 'Experiment' not found.")
  
  gene_col <- pick_gene_col(lit)
  
  # one ALL-like set per group (no UP/DOWN)
  split_tbl <- split(lit, lit[[group_col]])
  lit_sets <- lapply(names(split_tbl), function(grp) {
    vals <- split_tbl[[grp]][[gene_col]]
    ens  <- to_ensembl(vals)
    if (!length(ens)) return(NULL)
    
    # Experiment label: if grouping by Chemical, suffix _ALL; else keep the Experiment value as-is
    exp_lab <- if (group_col == "Chemical") paste0(grp, "_ALL") else as.character(grp)
    
    data.frame(
      Feature    = unique(ens),
      timepoint  = 24L,
      Experiment = exp_lab,
      stringsAsFactors = FALSE
    )
  })
  names(lit_sets) <- names(split_tbl)
  GList_lit <- Filter(Negate(is.null), lit_sets)
  if (!length(GList_lit)) stop("No gene sets produced after mapping to Ensembl.")
  
  # normalize ENSG versions if helper exists
  if (exists("normalize_glist")) GList_lit <- normalize_glist(GList_lit)
  
  # run enrichment against literature background
  ke_lit <- enrich_ke(GList_lit, "Literature", background = bg_lit)
} else {
  stop("File not found: ", LIT_FILE)
}


# ke_lit <- NULL
# LIT_FILE <- "outputs/DEG_tables/Literature_DEGs_Expanded.xlsx"
# 
# if (file.exists(LIT_FILE)) {
#   lit <- readr::read_csv(LIT_FILE, show_col_types = FALSE) |>
#          dplyr::filter(grepl("human|homo sapiens", Organism, ignore.case = TRUE))
#   lit$Direction <- tolower(lit$Direction)
#   lit_sets <- list()
#   splitC <- split(lit, lit$Chemical)
#   for (chem in names(splitC)) {
#     if ("UP"   %in% RUN_ON) lit_sets[[paste0(chem,"_UP")]]   <- unique(toupper(splitC[[chem]]$Gene[splitC[[chem]]$Direction=="up"]))
#     if ("DOWN" %in% RUN_ON) lit_sets[[paste0(chem,"_DOWN")]] <- unique(toupper(splitC[[chem]]$Gene[splitC[[chem]]$Direction=="down"]))
#     if ("ALL"  %in% RUN_ON) lit_sets[[paste0(chem,"_ALL")]]  <- unique(toupper(splitC[[chem]]$Gene))
#   }
#   # SYMBOL -> ENSEMBL for AOPfingerprintR
#   map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(unlist(lit_sets, use.names = FALSE)),
#                                keytype = "SYMBOL", columns = "ENSEMBL")
#   GList_lit <- lapply(names(lit_sets), function(nm) {
#     ensg <- unique(na.omit(map$ENSEMBL[match(lit_sets[[nm]], map$SYMBOL)]))
#     if (length(ensg)) data.frame(Feature = ensg, timepoint = 24L, Experiment = nm) else NULL
#   })
#   names(GList_lit) <- names(lit_sets)
#   GList_lit <- Filter(Negate(is.null), GList_lit)
#   GList_lit <- normalize_glist(GList_lit)  # <- normalize here too
# 
#   ke_lit <- enrich_ke(GList_lit, "Literature", background = bg_lit)
# }

# ---- Compare KEs across chemicals (heatmaps of -log10 q) --------------------
mk_matrix_from_ke <- function(ke_obj) {
  if (is.null(ke_obj)) return(NULL)

  # 1) Flatten if list of data.frames
  df <- if (is.data.frame(ke_obj)) {
    ke_obj
  } else {
    parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, ke_obj)
    if (!length(parts)) return(NULL)
    dplyr::bind_rows(parts, .id = "InputSet")
  }
  if (!nrow(df)) return(NULL)

  # 2) Find KE column
  nm <- names(df)
  ke_col <- if ("Ke" %in% nm) "Ke" else if ("TermID" %in% nm) "TermID" else nm[1]

  # 3) Keep only needed columns
  keep <- intersect(c(ke_col, "Experiment", "InputSet", "padj"), names(df))
  df <- df[, keep, drop = FALSE]

  # 3b) Ensure required columns exist with correct types
  if (!"Experiment" %in% names(df)) df$Experiment <- NA_character_
  if (!"InputSet"   %in% names(df)) df$InputSet   <- NA_character_
  if (!"padj"       %in% names(df)) df$padj       <- NA_real_

  # 4) Standardize and clean
  df <- df |>
    dplyr::mutate(
      KE   = .data[[ke_col]],
      EXP  = dplyr::coalesce(.data[["Experiment"]], .data[["InputSet"]], "set"),
      padj = suppressWarnings(as.numeric(.data[["padj"]]))
    ) |>
    dplyr::filter(!is.na(KE), nzchar(KE), !is.na(EXP), is.finite(padj))

  if (!nrow(df)) return(NULL)

  # 5) Collapse duplicates
  df_collapsed <- df |>
    dplyr::group_by(KE, EXP) |>
    dplyr::summarise(padj = min(padj), .groups = "drop") |>
    dplyr::mutate(score = -log10(pmax(padj, .Machine$double.xmin))) |>
    dplyr::select(KE, EXP, score)

  if (!nrow(df_collapsed)) return(NULL)

  # 6) Wide matrix
  wide <- tidyr::pivot_wider(
    df_collapsed,
    names_from  = EXP,
    values_from = score,
    values_fill = 0
  )

  if (!nrow(wide) || ncol(wide) < 2) return(NULL)

  wide |>
    tibble::column_to_rownames("KE") |>
    as.matrix()
}

if (!is.null(ke_gse)) {
  M <- mk_matrix_from_ke(ke_gse)
  if (!is.null(M)) {
    readr::write_csv(as.data.frame(M) |> tibble::rownames_to_column("KE"),
                     file.path(AOP_DIR, "KE_matrix_GSE.csv"))
    pheatmap::pheatmap(M, color = colorRampPalette(c("white","firebrick"))(101),
                       main = "KE enrichment (GSE DEGs)  -log10(q)",
                       filename = file.path(AOP_DIR,"KE_heatmap_GSE.png"),
                       width = 10, height = 8)
  }
}
if (!is.null(ke_lit)) {
  M <- mk_matrix_from_ke(ke_lit)
  if (!is.null(M)) {
    readr::write_csv(as.data.frame(M) |> tibble::rownames_to_column("KE"),
                     file.path(AOP_DIR, "KE_matrix_Literature.csv"))
    pheatmap::pheatmap(M, color = colorRampPalette(c("white","steelblue"))(101),
                       main = "KE enrichment (Literature) -log10(q)",
                       filename = file.path(AOP_DIR,"KE_heatmap_Literature.png"),
                       width = 10, height = 8)
  }
}


# ---- Visualize KEKE Network ----------------------------------
# simple slugify fallback
slugify <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "", x)
  tolower(x)
}

# Full-screen visNetwork saver
save_visnetwork <- function(widget, base, selfcontained = TRUE) {
  stopifnot(inherits(widget, "htmlwidget"))
  html_file <- paste0(base, ".html")
  
  # Make the widget fill the browser window
  widget$width  <- "100%"
  widget$height <- "100vh"
  widget$sizingPolicy <- htmlwidgets::sizingPolicy(
    defaultWidth  = "100%",
    defaultHeight = "100vh",
    browser.fill  = TRUE,   # fill external browser
    viewer.fill   = TRUE,   # fill RStudio viewer (as much as it can)
    padding       = 0
  )
  
  # Remove default margins and ensure 100% height from <html>/<body>
  widget <- htmlwidgets::prependContent(
    widget,
    htmltools::tags$style(htmltools::HTML("
      html, body { height:100%; margin:0; }
      .html-widget { height:100vh !important; width:100vw !important; }
    "))
  )
  
  # Save (optionally non-selfcontained for huge graphs)
  if (isTRUE(selfcontained)) {
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = TRUE)
  } else {
    libdir <- paste0(base, "_libs")
    dir.create(libdir, showWarnings = FALSE, recursive = TRUE)
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = FALSE, libdir = libdir)
  }
}


kes <- list(ke_gse, ke_lit)

for (i in seq_along(kes)) {
  ke <- kes[[i]]
  
  expts <- unique(na.omit(ke$Experiment))
  if (length(expts) == 0) {
    message("[", i, "] no Experiment values found.")
    next
  }
  
  for (ex in expts) {
    sub_out <- file.path(AOP_DIR, paste0("network_", slugify(ex)))
    dir.create(sub_out, showWarnings = FALSE, recursive = TRUE)

    ke_id = "TermID"
    numerical_variables = NULL
    pval_variable = "padj"
    gene_variable = "Genes"
    experiment = ex
    enlarge_ke_selection = T
    group_by = "ssbd"
    
    nodes_edges = suppressPackageStartupMessages(
      AOPfingerprintR::make_visNetwork(ke,
                                       experiment = experiment,
                                       enlarge_ke_selection = enlarge_ke_selection,
                                       ke_id, numerical_variables = numerical_variables,
                                       pval_variable,
                                       gene_variable,
                                       max_path_length = 40,
                                       n_AOs = 30,
                                       n_MIEs = 30)
    )
    
    write.csv(nodes_edges$nodes, file.path(sub_out, "nodes.csv"), row.names = FALSE)
    write.csv(nodes_edges$edges, file.path(sub_out, "edges.csv"), row.names = FALSE)
    
    vn = suppressPackageStartupMessages(
      AOPfingerprintR::plot_visNetwork(nodes = nodes_edges$nodes,
                                       edges = nodes_edges$edges,
                                       group_by = group_by,
                                       numerical_variables = numerical_variables)
    )
    
    save_visnetwork(vn, file.path(sub_out, "network"))
  }
}


# # 7) Lightweight run summary
# sink(file.path(AOP_DIR, sprintf("%s_summary.txt", chem_slug)))
# cat("Chemical: ", chem, "\nFile: ", infile, "\n\n", sep = "")
# cat("Sheets loaded: ", paste(sheets, collapse = ", "), "\n")
# if (is.data.frame(ke_enrichment)) {
#   cat("KE rows: ", nrow(ke_enrichment), "\n")
# }
# if (is.data.frame(aop_enrichment)) {
#   cat("AOP rows: ", nrow(aop_enrichment), "\n")
# }
# if (!is.null(res$detailed_results_only_enriched)) {
#   cat("Enriched rows: ", nrow(res$detailed_results_only_enriched), "\n")
#   cat("Experiments (enriched): ",
#       paste(unique(res$detailed_results_only_enriched$Experiment), collapse = ", "), "\n")
#   cat("AOPs (enriched): ", length(unique(res$detailed_results_only_enriched$AOP_name)), "\n")
# } else {
#   cat("No enriched AOPs/KEs at the selected thresholds.\n")
# }
# cat("\nSession info (short):\n")
# print(utils::sessionInfo()$R.version[c("version.string")])
# sink()
# 
# message("✓ Done: ", chem)





# ---- Map to AOs and draw an AO fingerprint ----------------------------------
make_fingerprint <- function(gene_GList, ke_df_or_list, label) {
  gene_GList <- lapply(gene_GList, function(df) {
    df <- df[, intersect(c("Feature","timepoint","Experiment"), names(df)), drop = FALSE]
    df
  })
  
  if (!length(gene_GList)) return(invisible(NULL))

  # --- background per label
  aop_bg <- if (tolower(label) == "literature") bg_lit else bg_gse

  # --- AOP enrichment
  aop_sig <- quiet_stdout(
    suppressPackageStartupMessages(enrich_KEs_AOPs(
      GList            = gene_GList,
      list_gene_sets   = human_ens_aop,
      only_significant = TRUE,
      pval_th          = 0.05,
      adj.method       = "fdr",
      merge_by         = "Aop",
      background       = aop_bg,
      numerical_properties = NULL,
      verbose          = FALSE
    ))
  )
  if (is.null(aop_sig)) {
    message("No significant AOPs for ", label); return(invisible(NULL))
  }

  # --- helper to bind list outputs
  .bind_df <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) return(x)
    parts <- Filter(function(d) is.data.frame(d) && nrow(d) > 0, x)
    if (!length(parts)) return(NULL)
    dplyr::bind_rows(parts, .id = "InputSet")
  }

  aop_df <- .bind_df(aop_sig)
  if (is.null(aop_df) || !nrow(aop_df)) {
    message("AOP table empty for ", label); return(invisible(NULL))
  }

  # --- KE table (needed by builder)
  ke_df <- .bind_df(ke_df_or_list)
  if (is.null(ke_df) || !nrow(ke_df)) {
    message("KE results empty; skip fingerprint for ", label); return(invisible(NULL))
  }

  # --- Build detailed fingerprint
  fp <- build_aop_for_aop_fingeprints(
    aop_enrichment_results = aop_df,
    ke_enrichment_results  = ke_df,
    min_aop_length         = 1,
    percentage_enriched_ke = 0.33
  )
  dr <- fp$detailed_results_only_enriched
  dr <- .bind_df(dr)
  if (is.null(dr) || !nrow(dr)) {
    message("No enriched AOP fingerprint rows for ", label); return(invisible(NULL))
  }

  # ------------------ STANDARDIZE / SAFEGUARD COLUMNS ------------------
  # AOP name
  if (!"AOP_name" %in% names(dr)) {
    if ("Aop" %in% names(dr))       dr <- dplyr::rename(dr, AOP_name = .data$Aop)
    else if ("AOP" %in% names(dr))  dr <- dplyr::rename(dr, AOP_name = .data$AOP)
    else                            dr$AOP_name <- paste0("AOP_", seq_len(nrow(dr)))
  }
  dr$AOP_name <- as.character(dr$AOP_name)

  # Experiment
  if (!"Experiment" %in% names(dr)) {
    cand <- intersect(c("InputSet","EXP","set","Group","Condition"), names(dr))
    dr$Experiment <- if (length(cand)) as.character(dr[[cand[1]]]) else "set"
  } else {
    dr$Experiment <- as.character(dr$Experiment)
  }

  # Grouping (SSbD_category) optional
  if (!"SSbD_category" %in% names(dr)) dr$SSbD_category <- "Uncategorized"
  dr$SSbD_category[is.na(dr$SSbD_category) | dr$SSbD_category==""] <- "Uncategorized"

  # Numeric size metric: prefer -log10(adjusted p); fallback to proportion enriched
  size_col <- NULL
  if ("padj" %in% names(dr)) {
    dr$neglog10_padj <- suppressWarnings(-log10(pmax(as.numeric(dr$padj), .Machine$double.xmin)))
    size_col <- "neglog10_padj"
  } else if ("p.adjust" %in% names(dr)) {
    dr$neglog10_padj <- suppressWarnings(-log10(pmax(as.numeric(dr$`p.adjust`), .Machine$double.xmin)))
    size_col <- "neglog10_padj"
  } else if ("pvalue" %in% names(dr)) {
    dr$neglog10_padj <- suppressWarnings(-log10(pmax(as.numeric(dr$pvalue), .Machine$double.xmin)))
    size_col <- "neglog10_padj"
  } else if ("proportion_enriched_KEs" %in% names(dr)) {
    dr$prop_enriched <- suppressWarnings(as.numeric(dr$proportion_enriched_KEs))
    size_col <- "prop_enriched"
  } else if ("percentage_enriched_ke" %in% names(dr)) {
    dr$prop_enriched <- suppressWarnings(as.numeric(dr$percentage_enriched_ke))/100
    size_col <- "prop_enriched"
  } else {
    # last resort: every point same small size
    dr$prop_enriched <- 0.2
    size_col <- "prop_enriched"
  }

  # Write standardized detailed table
  readr::write_csv(dr, file.path(AOP_DIR, paste0("AOP_fingerprint_", label, "_enriched.csv")))

  # ------------------ SAFE BUBBLE PLOT (NO PACKAGE PLOTTER) ------------------
  if (!all(c("AOP_name","Experiment", size_col) %in% names(dr))) {
    message("Required columns for plotting are missing; skipping plot for ", label)
    return(invisible(NULL))
  }

  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(dr, aes(x = Experiment, y = AOP_name)) +
    geom_point(aes(size = .data[[size_col]], fill = SSbD_category),
               shape = 21, alpha = 0.85, stroke = 0.2) +
    scale_size_continuous(name = ifelse(size_col == "neglog10_padj", "-log10(q)", "Proportion enriched"),
                          range = c(1.5, 10)) +
    guides(fill = guide_legend(title = "SSbD Category"), size = guide_legend(override.aes = list(shape = 21))) +
    labs(title = paste0("AOP fingerprint (", label, ")"),
         x = "Experiment / Input Set", y = "AOP") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")

  ggplot2::ggsave(file.path(AOP_DIR, paste0("AOP_fingerprint_", label, ".png")),
                  p, width = 10, height = 7, dpi = 600, bg = "white")
}

if (!is.null(ke_gse)) make_fingerprint(GList_gse, ke_gse, "GSE")
if (!is.null(ke_lit)) make_fingerprint(GList_lit, ke_lit, "Literature")

dir.create(file.path(AOP_DIR, "_provenance"), showWarnings = FALSE, recursive = TRUE)
readr::write_lines(bg_gse, file.path(AOP_DIR, "_provenance", "background_GSE_AS_IS.txt"))
readr::write_lines(bg_lit, file.path(AOP_DIR, "_provenance", "background_LIT_AS_IS.txt"))


cat("✓ Enrichment + comparisons complete. See:", normalizePath(AOP_DIR), "\n")

# ---- session info -------------------------------------------------------------
if (interactive()) {
    cat("\nSession info:\n")
    print(sessionInfo())
}

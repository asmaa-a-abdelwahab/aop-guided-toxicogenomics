#!/usr/bin/env Rscript
# ======================================================================
# Generic Functional Enrichment Pipeline
# - Case 1: DE-analysis DEGs (no TF expansion)
# - Case 2: Literature DEGs with TF expansion (ABC & AB; ALL/UP/DOWN; expanded/original)
#
# Outputs:
#   outputs/enrichment/{GO_BP,GO_MF,GO_CC,KEGG,Reactome}/...
#   outputs/enrichment/_summaries/*.csv (set sizes, overlaps, enrichment overlaps)
#   outputs/enrichment/all_enrich_objs.rds
#
# Assumes DE files under outputs/DE as produced by your DESeq2 script:
#   DEG_{CHEM}_vs_CTRL_all.csv
#   DEG_{CHEM}_{UP|DOWN|ALL}_ensg.txt (optional)
#
# Assumes TF expansion files under outputs/literature_tf_expansion (from DoRothEA script):
#   Chem_[ALL|UP|DOWN]_{ABC|AB}_{expanded|original}.csv
# ======================================================================

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "dplyr", "readr", "tidyr", "stringr", "tibble", "ggplot2",
  "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "ReactomePA",
  "enrichplot", "msigdbr"
))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr); library(tibble)
  library(ggplot2); library(grDevices)
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(enrichplot)
})

# ----------------------- CONFIG (edit here) ---------------------------
CFG <- list(
  # CASE 1 (DE-analysis) inputs
  DE_DIR        = DE_DIR,
  
  
  # NEW: turn DE-analysis on/off entirely
  RUN_DE        = TRUE,                     # <- set FALSE to skip DE-analysis block
  
  # CASE 2 (Literature + TF expansion) inputs
  TF_EXP_DIR    = file.path(OUT_DIR, "literature_tf_expansion", "Literature_DEGs"),
  TF_INCLUDE    = TRUE,                  # set FALSE to skip TF-expansion case
  TF_CONF       = c("ABC","AB"),        # which confidence bundles to include
  TF_PARTS      = c("ALL","UP","DOWN"), # which subsets to include
  TF_LIST_KIND  = c("expanded","original"), # which list kinds to enrich
  TF_CHEMS_FILTER = NULL, 
  
  # Optional limit/filter for DE contrasts (e.g., process only 5)
  CHEMS_FILTER  = NULL,                 # e.g., c("Ta4C3","Mo2Ti2C3","Nb4C3","Ti3C2Tx","LPS")
  MAX_CONTRASTS = Inf,
  
  
  # Enrichment parameters
  ENR_ROOT      = ENR_DIR,
  RUN_ON        = c("UP","DOWN","ALL"),
  Q_CUTOFF      = 0.05,
  MIN_GS        = 10,
  MAX_GS        = 500,
  TOP_N_PLOT    = 20,
  
  # Plotting
  PLOT_DPI      = 800,
  PLOT_WIDTH    = 20,
  PLOT_HEIGHT   = 12
)

# ----------------------- UTILITIES ------------------------------------
strip_v <- function(x) sub("\\.\\d+$", "", as.character(x))

ensure_dirs <- function(root) {
  dirs <- file.path(root, c("GO_BP","GO_MF","GO_CC","KEGG","Reactome","_summaries"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

save_both <- function(plot, basepath, width = CFG$PLOT_WIDTH, height = CFG$PLOT_HEIGHT, dpi = CFG$PLOT_DPI) {
  if (is.null(plot)) return(invisible(NULL))
  ggsave(paste0(basepath, ".png"), plot, width = width, height = height, dpi = dpi)
  if (capabilities("cairo")) {
    ggsave(paste0(basepath, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf)
  } else {
    ggsave(paste0(basepath, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, useDingbats = FALSE)
  }
}

wrap_descriptions <- function(x, width = 60) {
  if (!is.null(x) && !is.null(x@result) && "Description" %in% colnames(x@result)) {
    x@result$Description <- stringr::str_wrap(x@result$Description, width = width)
  }
  x
}

to_numeric_ratio <- function(x) {
  x <- as.character(x)
  vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
}

# SYMBOL -> ENSEMBL (unique, NA-safe)
map_symbol_to_ensembl <- function(sym) {
  if (!length(sym)) return(character(0))
  m <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = unique(na.omit(toupper(sym))),
                          keytype = "SYMBOL",
                          columns = "ENSEMBL")
  )
  unique(na.omit(sub("\\.\\d+$", "", m$ENSEMBL)))
}

plot_enrich_bar <- function(obj, top_n = CFG$TOP_N_PLOT, value = c("GeneRatio","Count"), showCategory = NULL) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj); if (!NROW(df)) return(NULL)
  if (!is.null(showCategory)) top_n <- showCategory
  value <- match.arg(value)
  if (value == "GeneRatio") {
    df$GeneRatioNum <- to_numeric_ratio(df$GeneRatio)
    xvar <- "GeneRatioNum"; ylab <- "GeneRatio"
  } else {
    xvar <- "Count";        ylab <- "Gene Count"
  }
  k <- min(top_n, nrow(df))
  df2 <- df |>
    arrange(p.adjust, desc(Count)) |>
    head(k)
  ggplot(df2, aes(x = reorder(Description, .data[[xvar]]), y = .data[[xvar]])) +
    geom_col() + coord_flip() +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 12)
}

safe_emapplot <- function(obj, showCategory = 30) {
  if (is.null(obj) || nrow(as.data.frame(obj)) <= 1) return(NULL)
  obj2 <- tryCatch(enrichplot::pairwise_termsim(obj), error = function(e) NULL)
  if (is.null(obj2) || is.null(obj2@termsim) || NROW(obj2@termsim) <= 1) return(NULL)
  tryCatch(enrichplot::emapplot(obj2, showCategory = showCategory, layout = "kk"), error = function(e) NULL)
}

safe_upsetplot <- function(obj, n = 15) {
  if (is.null(obj) || nrow(as.data.frame(obj)) <= 1) return(NULL)
  tryCatch(enrichplot::upsetplot(obj, n = n), error = function(e) NULL)
}

map_ensg_to_entrez <- function(ensg_vec) {
  if (!length(ensg_vec)) return(integer(0))
  ensg_vec <- unique(na.omit(strip_v(ensg_vec)))
  m <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = ensg_vec, keytype = "ENSEMBL",
                             columns = "ENTREZID")
  unique(na.omit(m$ENTREZID))
}

# ----------------------- ENRICHMENT CORE --------------------------------
run_enrichment_for_gene_set <- function(set_name, ensg, universe_entrez = NULL, fc_map = NULL) {
  message("\n---- Enriching set: ", set_name, "  (n=", length(ensg), ")")
  if (!length(ensg)) return(invisible(NULL))
  genes_entrez <- map_ensg_to_entrez(ensg)
  if (!length(genes_entrez)) { message("    (no mappable ENTREZ)"); return(invisible(NULL)) }
  
  root <- CFG$ENR_ROOT
  pAdj <- "BH"; qCut <- CFG$Q_CUTOFF
  MIN  <- CFG$MIN_GS; MAX <- CFG$MAX_GS
  
  # GO helper
  do_go <- function(ont) {
    ego <- suppressMessages(
      enrichGO(gene = genes_entrez, universe = universe_entrez,
               OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
               pAdjustMethod = pAdj, qvalueCutoff = qCut,
               minGSSize = MIN, maxGSSize = MAX, readable = TRUE)
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego <- wrap_descriptions(ego, 60)
      base <- file.path(root, paste0("GO_", ont), paste0("GO_", ont, "_", set_name))
      write.csv(as.data.frame(ego), paste0(base, ".csv"), row.names = FALSE)
      sig_top <- as.data.frame(ego) %>% filter(p.adjust < 0.05) %>% head(50)
      write.csv(sig_top, paste0(base, "_significant.csv"), row.names = FALSE)
      save_both(enrichplot::dotplot(ego, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ego)))) +
                  ggtitle(paste0("GO-", ont, ": ", set_name)),
                paste0(base, "_dotplot"))
      save_both(plot_enrich_bar(ego, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ego)))) +
                  ggtitle(paste0("GO-", ont, " barplot: ", set_name)),
                paste0(base, "_barplot"))
      save_both(safe_emapplot(ego, showCategory = 30), paste0(base, "_emapplot"))
      
      if (!is.null(fc_map)) {
        ensg_for_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = names(fc_map),
                                                 keytype = "ENSEMBL", columns = "ENTREZID")
        lfc_entrez <- setNames(fc_map[match(ensg_for_entrez$ENSEMBL, names(fc_map))],
                               ensg_for_entrez$ENTREZID)
        save_both(tryCatch(enrichplot::cnetplot(ego, showCategory = 10, foldChange = lfc_entrez),
                           error = function(e) NULL),
                  paste0(base, "_cnetplot"))
      }
      save_both(safe_upsetplot(ego, n = 15), paste0(base, "_upsetplot"))
    }
    invisible(ego)
  }
  
  ego_bp <- do_go("BP"); ego_mf <- do_go("MF"); ego_cc <- do_go("CC")
  
  # KEGG (with msigdbr fallback)
  ekegg <- NULL
  for (attempt in 1:3) {
    ekegg <- tryCatch(
      enrichKEGG(gene = genes_entrez, organism = "hsa",
                 pAdjustMethod = pAdj, qvalueCutoff = qCut,
                 minGSSize = MIN, maxGSSize = MAX,
                 universe = universe_entrez),
      error = function(e) NULL
    )
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) break
    Sys.sleep(2^attempt)
  }
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    # --- Robust MSigDB fallback: use KEGG only if available in this msigdbr version ---
    suppressPackageStartupMessages(require(msigdbr, quietly = TRUE))
    get_msig_kegg <- function() {
      avail <- tryCatch(msigdbr::msigdbr_collections(species = "Homo sapiens"),
                        error = function(e) NULL)
      if (is.null(avail) || !"subcategory" %in% names(avail)) return(NULL)
      # Normalize for matching
      avail$category    <- toupper(avail$category)
      avail$subcategory <- toupper(avail$subcategory)
      has_c2_kegg <- any(avail$category == "C2" & avail$subcategory %in% c("CP:KEGG","KEGG"))
      if (!has_c2_kegg) return(NULL)
      
      # Prefer CP:KEGG if present; else KEGG
      subcat <- if ("CP:KEGG" %in% avail$subcategory[avail$category == "C2"]) "CP:KEGG" else "KEGG"
      df <- tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = subcat),
        error = function(e) NULL
      )
      if (is.null(df) || !all(c("gs_name","entrez_gene") %in% names(df))) return(NULL)
      list(
        term2gene = setNames(df[, c("gs_name","entrez_gene")], c("term","gene")),
        term2name = setNames(unique(df[, c("gs_name","gs_name")]), c("term","name")),
        subcat = subcat
      )
    }
    
    msig <- get_msig_kegg()
    if (is.null(msig)) {
      message("  > KEGG: no significant enrichKEGG and no KEGG collection in msigdbr; skipping KEGG.")
    } else {
      ekegg <- tryCatch(
        enricher(gene = genes_entrez, universe = universe_entrez,
                 TERM2GENE = msig$term2gene, TERM2NAME = msig$term2name,
                 pAdjustMethod = pAdj, qvalueCutoff = qCut,
                 minGSSize = MIN, maxGSSize = MAX),
        error = function(e) NULL
      )
      if (!is.null(ekegg)) {
        message("  > KEGG via MSigDB (", msig$subcat, ") used as fallback.")
      }
    }
  }
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_r <- tryCatch(setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                        error = function(e) ekegg)
    ekegg_r <- wrap_descriptions(ekegg_r, 60)
    base <- file.path(CFG$ENR_ROOT, "KEGG", paste0("KEGG_", set_name))
    write.csv(as.data.frame(ekegg_r), paste0(base, ".csv"), row.names = FALSE)
    sig_top <- as.data.frame(ekegg_r) %>% filter(p.adjust < 0.05) %>% head(50)
    write.csv(sig_top, paste0(base, "_significant.csv"), row.names = FALSE)
    save_both(enrichplot::dotplot(ekegg_r, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ekegg_r)))),
              paste0(base, "_dotplot"))
    save_both(plot_enrich_bar(ekegg_r, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ekegg_r)))),
              paste0(base, "_barplot"))
    save_both(safe_emapplot(ekegg_r, showCategory = 30), paste0(base, "_emapplot"))
    if (!is.null(fc_map)) {
      ensg_for_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = names(fc_map),
                                               keytype = "ENSEMBL", columns = "ENTREZID")
      lfc_entrez <- setNames(fc_map[match(ensg_for_entrez$ENSEMBL, names(fc_map))],
                             ensg_for_entrez$ENTREZID)
      save_both(tryCatch(enrichplot::cnetplot(ekegg_r, showCategory = 10, foldChange = lfc_entrez),
                         error = function(e) NULL),
                paste0(base, "_cnetplot"))
    }
    save_both(safe_upsetplot(ekegg_r, n = 15), paste0(base, "_upsetplot"))
  }
  
  # Reactome
  ereact <- suppressMessages(
    ReactomePA::enrichPathway(gene = genes_entrez, organism = "human",
                              pAdjustMethod = pAdj, qvalueCutoff = qCut,
                              minGSSize = MIN, maxGSSize = MAX,
                              universe = universe_entrez, readable = TRUE)
  )
  if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
    ereact <- wrap_descriptions(ereact, 60)
    base <- file.path(CFG$ENR_ROOT, "Reactome", paste0("Reactome_", set_name))
    write.csv(as.data.frame(ereact), paste0(base, ".csv"), row.names = FALSE)
    sig_top <- as.data.frame(ereact) %>% filter(p.adjust < 0.05) %>% head(50)
    write.csv(sig_top, paste0(base, "_significant.csv"), row.names = FALSE)
    save_both(enrichplot::dotplot(ereact, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ereact)))),
              paste0(base, "_dotplot"))
    save_both(plot_enrich_bar(ereact, showCategory = min(CFG$TOP_N_PLOT, nrow(as.data.frame(ereact)))),
              paste0(base, "_barplot"))
    save_both(safe_emapplot(ereact, showCategory = 30), paste0(base, "_emapplot"))
    if (!is.null(fc_map)) {
      ensg_for_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = names(fc_map),
                                               keytype = "ENSEMBL", columns = "ENTREZID")
      lfc_entrez <- setNames(fc_map[match(ensg_for_entrez$ENSEMBL, names(fc_map))],
                             ensg_for_entrez$ENTREZID)
      save_both(tryCatch(enrichplot::cnetplot(ereact, showCategory = 10, foldChange = lfc_entrez),
                         error = function(e) NULL),
                paste0(base, "_cnetplot"))
    }
    save_both(safe_upsetplot(ereact, n = 15), paste0(base, "_upsetplot"))
  }
  
  # return objects (to collect if desired)
  list(GO_BP = ego_bp, GO_MF = ego_mf, GO_CC = ego_cc, KEGG = ekegg, Reactome = ereact)
}

# ----------------------- COMPARISON HELPERS ----------------------------
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

compare_gene_sets <- function(name_a, genes_a, name_b, genes_b) {
  tibble(
    setA = name_a,
    setB = name_b,
    nA = length(unique(genes_a)),
    nB = length(unique(genes_b)),
    n_intersect = length(intersect(genes_a, genes_b)),
    jaccard = jaccard(genes_a, genes_b)
  )
}

compare_enrichment_terms <- function(pathA, pathB, out_csv) {
  if (!file.exists(pathA) || !file.exists(pathB)) return(invisible(NULL))
  a <- tryCatch(readr::read_csv(pathA, show_col_types = FALSE), error = function(e) NULL)
  b <- tryCatch(readr::read_csv(pathB, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(a) || is.null(b) || !"ID" %in% names(a) || !"ID" %in% names(b)) return(invisible(NULL))
  df <- tibble(
    nA = nrow(a), nB = nrow(b),
    n_intersect = length(intersect(a$ID, b$ID)),
    jaccard = jaccard(a$ID, b$ID)
  )
  write_csv(df, out_csv)
  invisible(df)
}

# ----------------------- DISCOVER DE CONTRASTS -------------------------
discover_de_contrasts <- function() {
  all_csv <- list.files(CFG$DE_DIR, pattern = "^DEG_.*_vs_CTRL_all\\.csv$", full.names = TRUE)
  if (!length(all_csv)) return(character(0))
  get_name <- function(p) sub("^DEG_(.*)_vs_CTRL_all\\.csv$", "\\1", basename(p))
  setNames(all_csv, vapply(all_csv, get_name, character(1)))
}

# Preferred DE gene lists (UP/DOWN/ALL)
get_de_gene_list <- function(cn, mode, sig_df = NULL) {
  p <- file.path(CFG$DE_DIR, sprintf("DEG_%s_%s_ensg.txt", cn, mode))
  if (file.exists(p)) return(unique(na.omit(strip_v(readr::read_lines(p)))))
  if (is.null(sig_df)) return(character(0))
  if (!"ENSEMBL" %in% names(sig_df)) sig_df$ENSEMBL <- strip_v(rownames(sig_df))
  if (mode == "UP")   ens <- unique(sig_df$ENSEMBL[sig_df$log2FoldChange > 0])
  if (mode == "DOWN") ens <- unique(sig_df$ENSEMBL[sig_df$log2FoldChange < 0])
  if (mode == "ALL")  ens <- unique(sig_df$ENSEMBL)
  unique(na.omit(strip_v(ens)))
}

# ----------------------- MAIN RUNNER -----------------------------------
main <- function() {
  ensure_dirs(CFG$ENR_ROOT)
  summaries_dir <- file.path(CFG$ENR_ROOT, "_summaries")
  objs <- list()
  
  # ===== CASE 1: DE-analysis DEGs =====
  contrasts <- discover_de_contrasts()
  if (isTRUE(CFG$RUN_DE) && length(contrasts)) {
    if (!is.null(CFG$CHEMS_FILTER)) {
      contrasts <- contrasts[intersect(names(contrasts), CFG$CHEMS_FILTER)]
    }
    if (is.finite(CFG$MAX_CONTRASTS) && length(contrasts) > CFG$MAX_CONTRASTS) {
      contrasts <- contrasts[seq_len(CFG$MAX_CONTRASTS)]
    }
    
    for (cn in names(contrasts)) {
      message("\n====================================")
      message("DE-analysis enrichment: ", cn)
      all_df <- suppressMessages(readr::read_csv(contrasts[[cn]], show_col_types = FALSE))
      if (!"ENSEMBL" %in% names(all_df)) {
        if (!is.null(rownames(all_df))) all_df$ENSEMBL <- strip_v(rownames(all_df)) else
          stop("ENSEMBL column not found in: ", contrasts[[cn]])
      }
      universe <- map_ensg_to_entrez(all_df$ENSEMBL)
      
      fc_map <- NULL
      if ("log2FoldChange" %in% names(all_df)) {
        tmp <- all_df[, c("ENSEMBL","log2FoldChange")]
        tmp$ENSEMBL <- strip_v(tmp$ENSEMBL)
        fc_map <- setNames(tmp$log2FoldChange, tmp$ENSEMBL)
      }
      
      sig_csv <- file.path(CFG$DE_DIR, sprintf("DEG_%s_vs_CTRL_significant.csv", cn))
      sig_df <- if (file.exists(sig_csv)) suppressMessages(readr::read_csv(sig_csv, show_col_types = FALSE)) else NULL
      
      for (mode in CFG$RUN_ON) {
        ensg <- get_de_gene_list(cn, mode, sig_df)
        if (!length(ensg)) { message("  > ", mode, ": (no genes)"); next }
        set_name <- paste(cn, mode, sep = "_")
        objs[[set_name]] <- run_enrichment_for_gene_set(set_name, ensg, universe_entrez = universe, fc_map = fc_map)
      }
    }
  } else if (!isTRUE(CFG$RUN_DE)) {
    message("Skipping DE-analysis enrichment (CFG$RUN_DE=FALSE).")
  } else {
    message("No DE contrasts found under ", normalizePath(CFG$DE_DIR))
  }
  
  
  # ===== CASE 2: Literature + TF-expansion =====
  if (isTRUE(CFG$TF_INCLUDE) && dir.exists(CFG$TF_EXP_DIR)) {
    message("\n====================================")
    message("TF-expansion enrichment from: ", CFG$TF_EXP_DIR)
    
    # Accept both original and expanded; expanded has AB/ABC; filenames can be either
    #   <chem>_genes_original.csv
    #   <chem>_genes_expanded_AB.csv
    #   <chem>_genes_expanded_ABC.csv
    # (Also supports legacy "<chem>_<PART>_<CONF>_<which>.csv")
    tf_files <- list.files(CFG$TF_EXP_DIR,
                           pattern = "(original|expanded(_AB|_ABC)?|_(AB|ABC)_(expanded|original))\\.csv$",
                           full.names = TRUE,
                           ignore.case = TRUE)
    
    if (!length(tf_files)) {
      message("No TF files found in ", CFG$TF_EXP_DIR)
    } else {
      tf_tbl <- tibble(path = tf_files) %>%
        mutate(base = basename(path),
               # which: original vs expanded
               which = dplyr::case_when(
                 grepl("original\\.csv$", base, ignore.case = TRUE) ~ "original",
                 grepl("expanded", base, ignore.case = TRUE) ~ "expanded",
                 TRUE ~ NA_character_
               ),
               # conf: AB/ABC only for expanded (infer from either suffix or embedded tag)
               conf = dplyr::case_when(
                 grepl("_ABC(\\.|_)", base, ignore.case = TRUE) | grepl("_ABC\\.csv$", base, ignore.case = TRUE) ~ "ABC",
                 grepl("_AB(\\.|_)",  base, ignore.case = TRUE) | grepl("_AB\\.csv$",  base, ignore.case = TRUE)  ~ "AB",
                 TRUE ~ NA_character_
               ),
               # part: if your filenames carry _ALL/_UP/_DOWN; else default "ALL"
               part = dplyr::case_when(
                 grepl("_ALL_",  base, ignore.case = TRUE) ~ "ALL",
                 grepl("_UP_",   base, ignore.case = TRUE) ~ "UP",
                 grepl("_DOWN_", base, ignore.case = TRUE) ~ "DOWN",
                 TRUE ~ "ALL"
               ),
               # chem: first field after the TFEXP_ prefix
               chem = sub("_.*$", "", sub("^TFEXP_", "", base, ignore.case = TRUE)))
      
      # Keep only kinds requested
      tf_tbl <- tf_tbl %>%
        dplyr::filter(which %in% CFG$TF_LIST_KIND,
                      part  %in% CFG$TF_PARTS)
      
      # For expanded lists, honor TF_CONF; for originals, keep regardless of conf
      if (!is.null(CFG$TF_CONF)) {
        tf_tbl <- tf_tbl %>%
          dplyr::filter(which == "original" | (!is.na(conf) & conf %in% CFG$TF_CONF))
      }
      
      # Optional chemical filter
      if (!is.null(CFG$TF_CHEMS_FILTER)) {
        tf_tbl <- dplyr::filter(tf_tbl, chem %in% CFG$TF_CHEMS_FILTER)
      }
      
      message("TF files discovered: ", nrow(tf_tbl))
      if (nrow(tf_tbl)) message(paste(" -", tf_tbl$base, collapse = "\n"))
      
      # ---- Run enrichment for each TF list -------------------------------------
      for (i in seq_len(nrow(tf_tbl))) {
        row <- tf_tbl[i,]
        
        df <- suppressMessages(readr::read_csv(row$path, show_col_types = FALSE))
        
        # Accept ENSEMBL/Feature/Gene/SYMBOL; map to ENSEMBL if needed
        gene_col <- dplyr::case_when(
          "ENSEMBL" %in% names(df) ~ "ENSEMBL",
          "Feature" %in% names(df) ~ "Feature",
          "Gene"    %in% names(df) ~ "Gene",
          "SYMBOL"  %in% names(df) ~ "SYMBOL",
          TRUE ~ NA_character_
        )
        if (is.na(gene_col)) { message("Skipping (no gene column): ", row$base); next }
        
        vals <- unique(na.omit(df[[gene_col]]))
        if (gene_col %in% c("ENSEMBL","Feature")) {
          ensg <- unique(na.omit(sub("\\.\\d+$", "", as.character(vals))))
        } else {
          ensg <- map_symbol_to_ensembl(vals)
        }
        if (!length(ensg)) { message("No mappable ENSEMBL in: ", row$base); next }
        
        # Build a readable set name; originals get a fixed tag
        conf_tag <- ifelse(row$which == "expanded",
                           ifelse(is.na(row$conf), "EXP", row$conf),
                           "ORIG")
        set_name <- paste0("TFEXP_", row$chem, "_", row$part, "_", conf_tag, "_", toupper(row$which))
        
        # Universe left NULL (library defaults); adjust if you want a global universe
        objs[[set_name]] <- run_enrichment_for_gene_set(set_name, ensg,
                                                        universe_entrez = NULL, fc_map = NULL)
      }
      
      # ----- Comparisons ---------------------------------------------------------
      # A) AB vs ABC (same chem/part, expanded only)
      comp_rows <- list()
      for (chem in unique(tf_tbl$chem)) {
        for (part in unique(tf_tbl$part)) {
          p_ABC <- list.files(
            CFG$TF_EXP_DIR, full.names = TRUE, ignore.case = TRUE,
            pattern = paste0("^TFEXP_", chem, ".*", part, ".*(expanded_ABC|ABC_expanded).*\\.csv$")
          )
          p_AB <- list.files(
            CFG$TF_EXP_DIR, full.names = TRUE, ignore.case = TRUE,
            pattern = paste0("^TFEXP_", chem, ".*", part, ".*(expanded_AB|AB_expanded).*\\.csv$")
          )
          if (length(p_ABC) && length(p_AB)) {
            a <- suppressMessages(readr::read_csv(p_ABC[1], show_col_types = FALSE))
            b <- suppressMessages(readr::read_csv(p_AB[1],  show_col_types = FALSE))
            pick <- function(df) if ("ENSEMBL" %in% names(df)) "ENSEMBL" else if ("Feature" %in% names(df)) "Feature" else if ("Gene" %in% names(df)) "Gene" else "SYMBOL"
            gA <- unique(na.omit(if (pick(a) %in% c("ENSEMBL","Feature")) sub("\\.\\d+$","",a[[pick(a)]]) else map_symbol_to_ensembl(a[[pick(a)]])))
            gB <- unique(na.omit(if (pick(b) %in% c("ENSEMBL","Feature")) sub("\\.\\d+$","",b[[pick(b)]]) else map_symbol_to_ensembl(b[[pick(b)]])))
            comp_rows[[length(comp_rows)+1]] <- compare_gene_sets(
              paste0("TFEXP_", chem, "_", part, "_ABC_EXPANDED"), gA,
              paste0("TFEXP_", chem, "_", part, "_AB_EXPANDED"),  gB
            )
            bpA <- file.path(CFG$ENR_ROOT, "GO_BP", paste0("GO_BP_TFEXP_", chem, "_", part, "_ABC_EXPANDED.csv"))
            bpB <- file.path(CFG$ENR_ROOT, "GO_BP", paste0("GO_BP_TFEXP_", chem, "_", part, "_AB_EXPANDED.csv"))
            compare_enrichment_terms(bpA, bpB,
                                     file.path(CFG$ENR_ROOT, "_summaries",
                                               paste0("Overlap_GO_BP_TFEXP_", chem, "_", part, "_EXPANDED_AB_vs_ABC.csv")))
          }
        }
      }
      if (length(comp_rows)) {
        bind_rows(comp_rows) %>% write_csv(file.path(CFG$ENR_ROOT, "_summaries",
                                                     "TFEXP_AB_vs_ABC_gene_set_overlap.csv"))
      }
      
      # B) expanded vs original (same chem/part; conf = AB and ABC separately if present)
      comp_rows2 <- list()
      for (chem in unique(tf_tbl$chem)) {
        for (part in unique(tf_tbl$part)) {
          p_orig <- list.files(
            CFG$TF_EXP_DIR, full.names = TRUE, ignore.case = TRUE,
            pattern = paste0("^TFEXP_", chem, ".*(genes_)?original.*\\.csv$")
          )
          if (!length(p_orig)) next
          # compare vs AB
          p_exp_ab <- list.files(
            CFG$TF_EXP_DIR, full.names = TRUE, ignore.case = TRUE,
            pattern = paste0("^TFEXP_", chem, ".*", part, ".*(expanded_AB|AB_expanded).*\\.csv$")
          )
          # compare vs ABC
          p_exp_abc <- list.files(
            CFG$TF_EXP_DIR, full.names = TRUE, ignore.case = TRUE,
            pattern = paste0("^TFEXP_", chem, ".*", part, ".*(expanded_ABC|ABC_expanded).*\\.csv$")
          )
          
          load_ids <- function(p) {
            df <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
            pick <- if ("ENSEMBL" %in% names(df)) "ENSEMBL" else if ("Feature" %in% names(df)) "Feature" else if ("Gene" %in% names(df)) "Gene" else "SYMBOL"
            if (pick %in% c("ENSEMBL","Feature")) unique(na.omit(sub("\\.\\d+$","",df[[pick]])))
            else map_symbol_to_ensembl(df[[pick]])
          }
          
          gO <- load_ids(p_orig[1])
          
          if (length(p_exp_ab)) {
            gE <- load_ids(p_exp_ab[1])
            comp_rows2[[length(comp_rows2)+1]] <- compare_gene_sets(
              paste0("TFEXP_", chem, "_", part, "_ORIG"), gO,
              paste0("TFEXP_", chem, "_", part, "_AB_EXPANDED"), gE
            )
          }
          if (length(p_exp_abc)) {
            gE <- load_ids(p_exp_abc[1])
            comp_rows2[[length(comp_rows2)+1]] <- compare_gene_sets(
              paste0("TFEXP_", chem, "_", part, "_ORIG"), gO,
              paste0("TFEXP_", chem, "_", part, "_ABC_EXPANDED"), gE
            )
          }
        }
      }
      if (length(comp_rows2)) {
        bind_rows(comp_rows2) %>%
          write_csv(file.path(CFG$ENR_ROOT, "_summaries",
                              "TFEXP_EXPANDED_vs_ORIGINAL_gene_set_overlap.csv"))
      }
    }
  } else {
    message("Skipping TF-expansion (CFG$TF_INCLUDE=FALSE or dir missing).")
  }
  
  
  # Save objects handle (optional exploration later)
  saveRDS(objs, file = file.path(CFG$ENR_ROOT, "all_enrich_objs.rds"))
  cat("\n✓ Finished. Enrichment saved under: ", normalizePath(CFG$ENR_ROOT), "\n")
}

# --- helpers to merge *_significant.csv per collection -----------------
.parse_experiment_from_file <- function(path) {
  base <- basename(path)
  nm <- sub("\\.csv$", "", base)
  nm <- sub("^GO_(BP|MF|CC)_", "", nm)   # strip GO prefix
  nm <- sub("^KEGG_", "", nm)            # strip KEGG prefix
  nm <- sub("^Reactome_", "", nm)        # strip Reactome prefix
  nm <- sub("_significant$", "", nm)     # strip suffix
  nm
}

merge_significant_per_collection <- function(root, out_dir,
                                             collections = c("GO_BP","GO_MF","GO_CC","KEGG","Reactome")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (coll in collections) {
    coll_dir <- file.path(root, coll)
    files <- list.files(coll_dir, pattern = "_significant\\.csv$", full.names = TRUE)
    if (!length(files)) { message("No significant files found in ", coll_dir); next }
    
    parts <- lapply(files, function(f) {
      df <- tryCatch(suppressMessages(readr::read_csv(f, show_col_types = FALSE)), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
      df$Experiment <- .parse_experiment_from_file(f)
      df
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) { message("All significant tables empty in ", coll_dir); next }
    
    merged <- dplyr::bind_rows(parts)
    out_file <- file.path(out_dir, paste0("merged_", coll, "_significant.csv"))
    readr::write_csv(merged, out_file)
    message("Merged ", nrow(merged), " rows -> ", normalizePath(out_file))
  }
}

ensure_dirs(CFG$ENR_ROOT)
main()

merge_significant_per_collection(
  root    = CFG$ENR_ROOT,
  out_dir = file.path(CFG$ENR_ROOT, "_summaries")
)

# ===================== end of script =====================


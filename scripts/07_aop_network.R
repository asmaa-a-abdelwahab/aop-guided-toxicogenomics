#!/usr/bin/env Rscript
# AOP fingerprint -> AOP-Wiki network visualizer.
# Produces a publication PNG, an interactive HTML file, and subgraph tables.

source(file.path("scripts", "_config.R"))
assert_packages(c(
  "optparse", "readr", "readxl", "dplyr", "tidyr", "stringr", "tibble",
  "igraph", "tidygraph", "ggraph", "visNetwork", "scales", "htmlwidgets"
))

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(visNetwork)
  library(scales)
})

option_list <- list(
  make_option("--enrichment", type = "character",
              default = file.path(AOP_DIR, "AOP_fingerprint_Literature_enriched.csv"),
              help = "AOP fingerprint enrichment CSV [default: %default]"),
  make_option("--ke-metadata", type = "character",
              default = file.path(ROOT, "data", "aop_wiki", "aop_ke_ec.tsv"),
              help = "AOP-Wiki KE metadata TSV [default: %default]"),
  make_option("--relationships", type = "character",
              default = file.path(ROOT, "data", "aop_wiki", "aop_ke_ker.tsv"),
              help = "AOP-Wiki KE relationship TSV [default: %default]"),
  make_option("--roles", type = "character",
              default = file.path(ROOT, "data", "aop_wiki", "aop_ke_mie_ao.tsv"),
              help = "AOP-Wiki KE role TSV [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = file.path(OUT_DIR, "AOP_network"),
              help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

in_enrich <- opt$enrichment
in_ke_ec  <- opt[["ke-metadata"]]
in_ker    <- opt$relationships
in_mieao  <- opt$roles
out_dir   <- opt$outdir

assert_files(
  c(in_enrich, in_ke_ec, in_ker, in_mieao),
  "AOP network input files"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SAFE_MODE <- TRUE
MAX_NODES <- 500
MAX_EDGES <- 2000

alpha <- 0.05
min_neglog10 <- -log10(alpha)
label_enriched_only <- TRUE
top_degree_to_label <- 20

png_width  <- 3800
png_height <- 2600
png_res    <- 350

## --- Helpers -----------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

coalesce_cols <- function(df, candidates, new_name) {
  x <- NULL
  for (c in candidates) if (c %in% names(df)) { x <- df[[c]]; break }
  if (is.null(x)) return(df %>% mutate(!!new_name := NA_character_))
  df %>% mutate(!!new_name := x)
}

num_coalesce <- function(df, candidates, new_name) {
  x <- NULL
  for (c in candidates) if (c %in% names(df)) { x <- suppressWarnings(as.numeric(df[[c]])); break }
  if (is.null(x)) return(df %>% mutate(!!new_name := NA_real_))
  df %>% mutate(!!new_name := x)
}

read_aopwiki_tsv <- function(path, fallback_names) {
  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  first_field <- if (length(first_line)) {
    tolower(trimws(strsplit(first_line, "\t", fixed = TRUE)[[1]][[1]]))
  } else {
    ""
  }
  has_header <- first_field %in% c("aop", "aop_id", "aop.id", "aop_wiki_id")

  suppressMessages(readr::read_tsv(
    path,
    col_names = if (has_header) TRUE else fallback_names,
    show_col_types = FALSE,
    name_repair = "minimal",
    trim_ws = TRUE
  ))
}

## --- Load: enrichment results -------------------------------------------------
enrich_raw <- readr::read_csv(in_enrich, show_col_types = FALSE)

# Try to detect key columns with multiple aliases
enrich <- enrich_raw %>%
  # KE identifiers/names
  coalesce_cols(c("Ke_description","KE","KeyEvent","Key_Event","KE_Name","KE_name","key_event"), "KE_name") %>%
  coalesce_cols(c("Ke","KE_ID","KEID","KeyEventID","KE_Id","KE.id","KE_Wiki_ID","KE_AOPWiki_ID"), "KE_id") %>%
  # AOP / Pathway identifiers/names
  coalesce_cols(c("AOP_name","a.name","AOP","AOP_ID","AOP.id","AOPID","AOP_Title","Pathway","Pathway_Name"), "AOP_label") %>%
  coalesce_cols(c("TermID","AOP_Wiki_ID","AOP_ID_num","AOP_Number","AOP_No","AOP.ID"), "AOP_id") %>%
  # scores / significance
  num_coalesce(c("padj","p.adjust","adj.p","FDR"), "padj") %>%
  num_coalesce(c("neglog10_padj","score","NES","-log10(padj)","minuslog10padj"), "neglog10_padj") %>%
  mutate(
    neglog10_padj = ifelse(is.na(neglog10_padj) & !is.na(padj), -log10(padj), neglog10_padj),
    is_signif = ifelse(!is.na(neglog10_padj), neglog10_padj >= min_neglog10,
                       ifelse(!is.na(padj), padj <= alpha, FALSE))
  )

# Collect enriched KEs and AOPs
enriched_KEs  <- enrich %>%
  filter(is_signif) %>%
  transmute(KE_id, KE_name, neglog10_padj, padj, AOP_id, AOP_label) %>%
  distinct()

enriched_AOPs <- enrich %>%
  filter(is_signif) %>%
  transmute(AOP_id, AOP_label) %>%
  distinct()

message(">>> Enriched KEs n = ", nrow(enriched_KEs))
message(">>> Enriched AOPs n = ", nrow(enriched_AOPs))

## --- Load: AOP-Wiki TSVs -----------------------------------------------------
# Expecting:
# - ke_ec: at least KE_id, KE_name (or similar)
# - ke_ker: columns like AOP_id, KE_up_id (source), KE_down_id (target), KER_id, evidence, weight_of_evidence
# - mie_ao: columns mapping KE_id to role (MIE/KE/AO) and AOP_id/AOP_title
ke_ec <- read_aopwiki_tsv(
  in_ke_ec,
  c(
    "AOP_id", "KE_id", "Effect",
    "Ontology_1", "Ontology_id_1", "KE_name_primary",
    "Ontology_2", "Ontology_id_2", "KE_name"
  )
)
ke_ker <- read_aopwiki_tsv(
  in_ker,
  c("AOP_id", "from_id", "to_id", "KER_id", "Adjacency", "Evidence", "Quantitative")
)
mieao <- read_aopwiki_tsv(
  in_mieao,
  c("AOP_id", "KE_id", "Role", "KE_name")
)
if (!"KE_name_primary" %in% names(ke_ec)) ke_ec$KE_name_primary <- NA_character_

# Normalize key columns
ke_ec <- ke_ec %>%
  coalesce_cols(c("KE_id","KE.ID","KEID","KE_wiki_id","Event_ID"), "KE_id") %>%
  coalesce_cols(c("KE_name","Event","Event_Name","KE Title","KE_title"), "KE_name") %>%
  mutate(
    KE_id = as.character(KE_id),
    KE_name = dplyr::coalesce(
      as.character(KE_name),
      as.character(KE_name_primary),
      KE_id
    )
  ) %>%
  distinct(KE_id, .keep_all = TRUE)

ke_ker <- ke_ker %>%
  coalesce_cols(c("AOP_id","AOP.ID","AOP","AOP_wiki_id"), "AOP_id") %>%
  coalesce_cols(c("AOP_title","AOP_Title","Pathway","AOP_label"), "AOP_label") %>%
  coalesce_cols(c("from_id","KE_up_id","KE_upstream_id","from","source","KE_source"), "from_id") %>%
  coalesce_cols(c("to_id","KE_down_id","KE_downstream_id","to","target","KE_target"), "to_id") %>%
  coalesce_cols(c("KER_id","KER.ID","Relationship_ID","Rel_ID"), "KER_id") %>%
  mutate(AOP_id = as.character(AOP_id),
         from_id = as.character(from_id),
         to_id   = as.character(to_id))

mieao <- mieao %>%
  coalesce_cols(c("AOP_id","AOP.ID","AOP_wiki_id"), "AOP_id") %>%
  coalesce_cols(c("AOP_title","AOP_Title","AOP_label"), "AOP_label") %>%
  coalesce_cols(c("KE_id","KE.ID","Event_ID"), "KE_id") %>%
  coalesce_cols(c("KE_name","Event","Event_Name","KE Title","KE_title"), "KE_name") %>%
  coalesce_cols(c("Role","KE_role","Node_Type"), "Role") %>%
  mutate(AOP_id = as.character(AOP_id),
         KE_id  = as.character(KE_id))


# ---------------- PATCH: normalize + map names to IDs + diagnostics ----------------

# 0) Quick peek at available columns
message(">>> Enrichment cols: ", paste(names(enrich_raw), collapse=", "))
message(">>> ke_ec cols: ", paste(names(ke_ec), collapse=", "))
message(">>> ke_ker cols: ", paste(names(ke_ker), collapse=", "))
message(">>> mieao cols: ", paste(names(mieao), collapse=", "))

# 1) Clean whitespace / case on keys to improve matching
trimlower <- function(x) { tolower(trimws(as.character(x))) }

ke_ec <- ke_ec %>%
  mutate(
    KE_id   = as.character(KE_id),
    KE_name = ifelse(is.na(KE_name), KE_id, KE_name),
    KE_name_norm = trimlower(KE_name)
  )

mieao <- mieao %>%
  mutate(
    AOP_id     = as.character(AOP_id),
    AOP_label  = ifelse(is.na(AOP_label), AOP_id, AOP_label),
    AOP_label_norm = trimlower(AOP_label)
  )

# 2) Normalize enrichment names/labels and map to IDs when missing
enrich <- enrich %>%
  mutate(
    KE_name_norm   = ifelse(is.na(KE_name), NA_character_, trimlower(KE_name)),
    AOP_label_norm = ifelse(is.na(AOP_label), NA_character_, trimlower(AOP_label))
  )

# 2a) If KE_id missing but KE_name present, map via ke_ec
if (all(is.na(enrich$KE_id)) && any(!is.na(enrich$KE_name_norm))) {
  message(">>> Mapping KE_name → KE_id from ke_ec (KE_id missing in enrichment).")
  enrich <- enrich %>%
    left_join(ke_ec %>% select(KE_id, KE_name_norm), by = "KE_name_norm") %>%
    mutate(KE_id = coalesce(KE_id.x, KE_id.y)) %>%
    select(-KE_id.x, -KE_id.y)
}

# 2b) If AOP_id missing but label present, map via mieao
if (all(is.na(enrich$AOP_id)) && any(!is.na(enrich$AOP_label_norm))) {
  message(">>> Mapping AOP_label → AOP_id from mieao (AOP_id missing in enrichment).")
  enrich <- enrich %>%
    left_join(mieao %>% select(AOP_id, AOP_label_norm) %>% distinct(), by = "AOP_label_norm") %>%
    mutate(AOP_id = coalesce(AOP_id.x, AOP_id.y)) %>%
    select(-AOP_id.x, -AOP_id.y)
}

# 3) Recompute significance flags (in case numeric coercion created NAs)
enrich <- enrich %>%
  mutate(
    neglog10_padj = ifelse(is.na(neglog10_padj) & !is.na(padj), -log10(padj), neglog10_padj),
    is_signif = ifelse(!is.na(neglog10_padj), neglog10_padj >= min_neglog10,
                       ifelse(!is.na(padj), padj <= alpha, FALSE))
  )

# 4) If still nothing significant, pick top N as a fallback so you always get a graph
if (!any(enrich$is_signif, na.rm = TRUE)) {
  message(">>> No rows pass padj threshold; falling back to top 25 by −log10(padj).")
  enrich <- enrich %>%
    mutate(fallback_rank = rank(desc(coalesce(neglog10_padj, -log10(coalesce(padj, 1)))))) %>%
    mutate(is_signif = fallback_rank <= 25)
}

# 5) Rebuild enriched_KEs / enriched_AOPs with the (now) mapped IDs
enriched_KEs  <- enrich %>%
  filter(is_signif) %>%
  transmute(KE_id = as.character(KE_id),
            KE_name, neglog10_padj, padj,
            AOP_id = as.character(AOP_id), AOP_label) %>%
  distinct()

enriched_AOPs <- enrich %>%
  filter(is_signif) %>%
  transmute(AOP_id = as.character(AOP_id), AOP_label) %>%
  distinct()

message(">>> After mapping: enriched KEs n = ", nrow(enriched_KEs),
        " | unique KE_ids = ", length(unique(na.omit(enriched_KEs$KE_id))))
message(">>> After mapping: enriched AOPs n = ", nrow(enriched_AOPs),
        " | unique AOP_ids = ", length(unique(na.omit(enriched_AOPs$AOP_id))))

# 6) If we still have NA KE_id but a name, try a second pass strict join KE_name_norm
if (any(is.na(enriched_KEs$KE_id)) && any(!is.na(enriched_KEs$KE_name))) {
  message(">>> Second-pass KE_name → KE_id mapping on enriched rows.")
  enriched_KEs <- enriched_KEs %>%
    mutate(KE_name_norm = trimlower(KE_name)) %>%
    left_join(ke_ec %>% select(KE_id, KE_name_norm), by = "KE_name_norm") %>%
    mutate(KE_id = coalesce(KE_id.x, KE_id.y)) %>%
    select(-KE_id.x, -KE_id.y, -KE_name_norm)
}

# 7) Final sanity logs before building ker_sub
message(">>> Non-NA enriched KE_ids: ", sum(!is.na(enriched_KEs$KE_id)))
message(">>> Non-NA enriched AOP_ids:", sum(!is.na(enriched_AOPs$AOP_id)))



## --- Focus: only AOPs touched by enrichment ----------------------------------
touched_aops <- dplyr::bind_rows(
  enriched_AOPs %>% filter(!is.na(AOP_id)) %>% transmute(AOP_id),
  mieao %>% filter(KE_id %in% enriched_KEs$KE_id) %>% transmute(AOP_id)
) %>% distinct() %>% pull(AOP_id)

message(">>> Touched AOP IDs: ", length(touched_aops),
        " | matching KER AOP IDs: ",
        length(intersect(touched_aops, unique(ke_ker$AOP_id))))

if (length(touched_aops) == 0) {
  warning("No AOP IDs detected from results; falling back to all KERs that include enriched KE IDs.")
}

# --- smart KER subset selection ---
ker_sub <- NULL

if (length(touched_aops)) {
  ker_sub <- ke_ker %>% filter(AOP_id %in% touched_aops)
}

# If still empty, fall back to KE_id membership in edges
if (is.null(ker_sub) || nrow(ker_sub) == 0) {
  message(">>> No KERs via AOP IDs; falling back to KERs that touch enriched KE IDs.")
  valid_ke_ids <- unique(na.omit(enriched_KEs$KE_id))
  ker_sub <- ke_ker %>% filter(from_id %in% valid_ke_ids | to_id %in% valid_ke_ids)
}

# Last resort: try matching by KE names (rare but helpful if IDs differ in sources)
if (nrow(ker_sub) == 0) {
  message(">>> Still empty; attempting KE name-based match via ke_ec.")
  # Map ker endpoints to names then match any enriched KE_name
  ker_named <- ke_ker %>%
    left_join(ke_ec %>% select(KE_id, KE_name), by = c("from_id" = "KE_id")) %>%
    rename(from_name = KE_name) %>%
    left_join(ke_ec %>% select(KE_id, KE_name), by = c("to_id" = "KE_id")) %>%
    rename(to_name = KE_name)
  
  ker_sub <- ker_named %>%
    filter(from_name %in% enriched_KEs$KE_name | to_name %in% enriched_KEs$KE_name)
}

# Keep only rows with valid endpoints
ker_sub <- ker_sub %>% filter(!is.na(from_id), !is.na(to_id))
message(">>> ker_sub rows: ", nrow(ker_sub))

if (nrow(ker_sub) == 0) {
  warning("No matching KERs after all mapping strategies. Check that KE IDs/titles in your enrichment file correspond to the AOP-Wiki extracts.")
  # Still write empty CSVs to keep pipeline stable, but you can stop here if preferred:
  readr::write_csv(tibble(), file.path(out_dir, "subgraph_nodes.csv"))
  readr::write_csv(tibble(), file.path(out_dir, "subgraph_edges.csv"))
  quit(save = "no", status = 0)
}


# Node table = all KEs participating in ker_sub
node_ids <- sort(unique(c(ker_sub$from_id, ker_sub$to_id)))
nodes <- tibble(KE_id = node_ids) %>%
  left_join(ke_ec %>% select(KE_id, KE_name), by = "KE_id") %>%
  left_join(
    mieao %>%
      group_by(KE_id) %>%
      summarize(KE_name_role = dplyr::first(
                  KE_name[!is.na(KE_name) & nzchar(KE_name)],
                  default = NA_character_
                ),
                Role = paste(sort(unique(na.omit(Role))), collapse = "; "),
                AOP_ids = paste(sort(unique(na.omit(AOP_id))), collapse = "; "),
                AOP_labels = paste(sort(unique(na.omit(AOP_label))), collapse = "; "),
                .groups = "drop"),
    by = "KE_id"
  ) %>%
  mutate(
    KE_name = dplyr::coalesce(KE_name_role, KE_name, KE_id),
    Role = ifelse(is.na(Role), "KE", Role)
  ) %>%
  select(-KE_name_role)

# Edge table for graph
edges <- ker_sub %>%
  transmute(from = from_id, to = to_id,
            AOP_id, AOP_label,
            KER_id = ifelse(is.na(KER_id), paste0(from,"→",to), KER_id))

# Add enrichment flags to nodes
nodes <- nodes %>%
  left_join(enriched_KEs %>% select(KE_id, neglog10_padj, padj) %>%
              distinct(), by = "KE_id") %>%
  mutate(
    is_enriched = !is.na(neglog10_padj) & neglog10_padj >= min_neglog10,
    enrich_score = ifelse(is.na(neglog10_padj), 0, neglog10_padj),
    node_label = dplyr::coalesce(KE_name, KE_id)
  )

# Enrichment flags for edges via AOP membership (for coloring subgraphs)
edges <- edges %>%
  mutate(is_in_enriched_aop = ifelse(!is.na(AOP_id) & AOP_id %in% enriched_AOPs$AOP_id, TRUE, FALSE))

## --- Build graph --------------------------------------------------------------
g_tbl <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

# Centrality for optional labeling
g_tbl <- g_tbl %>%
  activate(nodes) %>%
  mutate(deg = centrality_degree(mode = "all"),
         btwn = centrality_betweenness(directed = TRUE, normalized = TRUE))

# --- Decide labels for PNG (safe with tbl_graph) ---
if (label_enriched_only) {
  # Get node data as a tibble to ensure a plain logical vector
  node_tbl <- as_tibble(g_tbl, active = "nodes")
  
  # Make sure it's logical with no NAs
  is_enriched_vec <- dplyr::coalesce(as.logical(node_tbl$is_enriched), FALSE)
  
  lab_ids <- which(is_enriched_vec)
  
  # also allow top-degree nodes (context)
  extra <- order(node_tbl$deg, decreasing = TRUE)[seq_len(min(top_degree_to_label, nrow(node_tbl)))]
  lab_ids <- sort(unique(c(lab_ids, extra)))
  
  g_tbl <- g_tbl %>%
    activate(nodes) %>%
    mutate(show_label = row_number() %in% lab_ids)
} else {
  g_tbl <- g_tbl %>%
    activate(nodes) %>%
    mutate(show_label = TRUE)
}

# Keep the largest weakly connected component to prevent layout blowups.
g_tbl <- tidygraph::convert(
  g_tbl,
  tidygraph::to_largest_component,
  type = "weak"
)

## --- VIS 1: Publication PNG (ggraph) -----------------------------------------
png_file <- file.path(out_dir, "AOP_enrichment_network.png")
message(">>> Writing PNG: ", png_file)

# Node aesthetics
node_fill <- ifelse(g_tbl$is_enriched, "enriched", "not_enriched")
node_shape <- dplyr::case_when(
  str_detect(g_tbl$Role %||% "", "AO")  ~ "triangle",
  str_detect(g_tbl$Role %||% "", "MIE") ~ "square",
  TRUE                                  ~ "circle"
)

# Build plot
set.seed(42)
# before building `p`, ensure node groups exist as factors
g_tbl <- g_tbl %>%
  activate(nodes) %>%
  mutate(
    fill_group  = factor(ifelse(is_enriched, "Enriched KE", "Other KE"),
                         levels = c("Enriched KE","Other KE")),
    shape_group = factor(dplyr::case_when(
      grepl("AO",  Role %||% "")  ~ "AO",
      grepl("MIE", Role %||% "")  ~ "MIE",
      TRUE                        ~ "KE"
    ), levels = c("KE","MIE","AO"))
  )

# rebuild plot using those factor columns
p <- ggraph(g_tbl, layout = "sugiyama") +
  geom_edge_link(aes(alpha = ifelse(is_in_enriched_aop, 0.45, 0.15)),
                 arrow = arrow(length = unit(3.5, "pt"), type = "closed"),
                 end_cap = circle(2.2, 'pt'), start_cap = circle(2.2, 'pt')) +
  geom_node_point(aes(size = 1 + sqrt(deg),
                      fill = fill_group,
                      shape = shape_group,
                      alpha = ifelse(is_enriched, 0.95, 0.75)),
                  stroke = 0.3, color = "grey20") +
  scale_fill_manual(values = c("Enriched KE"="#d62728", "Other KE"="#9ecae1")) +
  scale_shape_manual(values = c("KE"=21, "MIE"=22, "AO"=24)) +
  guides(alpha = "none") +
  theme_graph(base_family = "sans") +
  labs(title = "AOP-Wiki Network with Enriched KEs / Pathways",
       subtitle = "Red = enriched KE | Squares = MIE | Triangles = AO | Circles = KE",
       fill = "KE status", shape = "Node type", size = "Degree (√)") +
  theme(legend.position = "right")

# labels (no x/y needed; ggraph supplies them)
p <- p + ggraph::geom_node_text(
  data = function(d) dplyr::filter(d, show_label),
  mapping = aes(label = node_label),
  repel = TRUE, size = 2.8, fontface = 2
)

png(png_file, width = png_width, height = png_height, res = png_res)
print(p)
invisible(dev.off())

## --- VIS 2: Interactive HTML (visNetwork) ------------------------------------
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub("\"","&quot;",x, fixed = TRUE)
  x <- gsub("'", "&#39;",  x, fixed = TRUE)
  x
}

html_file <- file.path(out_dir, "AOP_enrichment_network.html")
message(">>> Writing HTML: ", html_file)

# Build data frames (same as you have)
vn_nodes <- as_tibble(g_tbl) %>%
  transmute(
    id    = KE_id,
    label = node_label,
    title = paste0(
      "<b>", html_escape(node_label), "</b>",
      ifelse(!is.na(Role), paste0("<br/>Role: ", html_escape(Role)), ""),
      ifelse(!is.na(AOP_labels) & AOP_labels != "",
             paste0("<br/>AOPs: ", html_escape(AOP_labels)), ""),
      ifelse(!is.na(padj), paste0("<br/>padj: ", signif(padj, 3)), ""),
      ifelse(!is.na(neglog10_padj), paste0("<br/>−log10(padj): ", round(neglog10_padj, 3)), "")
    ),
    shape = dplyr::case_when(
      grepl("AO",  Role %||% "") ~ "triangle",
      grepl("MIE", Role %||% "") ~ "square",
      TRUE                       ~ "dot"
    ),
    color = ifelse(is_enriched, "#d62728", "#6aaed6"),
    size  = pmax(8, pmin(30, 8 + 4*sqrt(deg))),
    is_enriched = is_enriched,
    deg = deg
  )

vn_edges <- as_tibble(activate(g_tbl, edges)) %>%
  transmute(from = from, to = to,
            title = paste0("AOP: ", AOP_label %||% AOP_id),
            arrows = "to",
            color = ifelse(is_in_enriched_aop, "#9E9E9E", "#D0D0D0"),
            width = ifelse(is_in_enriched_aop, 2.2, 1))

# Optional: cap size for interactive rendering (keep all enriched + top-degree context)
if (SAFE_MODE) {
  keep_ids <- vn_nodes %>%
    arrange(desc(is_enriched), desc(deg)) %>%
    slice_head(n = MAX_NODES) %>%
    pull(id)
  
  vn_nodes <- vn_nodes %>% filter(id %in% keep_ids)
  vn_edges <- vn_edges %>% filter(from %in% keep_ids & to %in% keep_ids)
  
  if (nrow(vn_edges) > MAX_EDGES) {
    vn_edges <- vn_edges %>% slice_sample(n = MAX_EDGES)
  }
}

html_file <- file.path(out_dir, "AOP_enrichment_network.html")
lib_dir   <- file.path(out_dir, "libs")
dir.create(lib_dir, showWarnings = FALSE)

# Build widget
vw <- visNetwork(vn_nodes, vn_edges, main = "AOP-Wiki Network (interactive)") %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE, solver = "forceAtlas2Based") %>%
  visLegend(useGroups = FALSE, addNodes = data.frame(
    label = c("Enriched KE","Non-enriched KE","MIE","AO"),
    shape = c("dot","dot","square","triangle"),
    color = c("#d62728","#6aaed6","#6aaed6","#6aaed6"),
    title = c("−log10(padj) ≥ 1.3", "Other KEs","MIE","AO")
  ))

# Save non-selfcontained to keep file small; wrap in tryCatch to avoid session crash
tryCatch(
  {
    htmlwidgets::saveWidget(vw, file = html_file, selfcontained = FALSE, libdir = lib_dir)
    message(">>> HTML saved (non-selfcontained).")
  },
  error = function(e) {
    message(">>> saveWidget failed: ", conditionMessage(e))
    # As a last resort, save a minimal self-contained (may be large)
    try(htmlwidgets::saveWidget(vw, file = html_file, selfcontained = TRUE), silent = TRUE)
  }
)


## --- Save subgraph tables for record -----------------------------------------
readr::write_csv(as_tibble(activate(g_tbl, nodes)), file.path(out_dir, "subgraph_nodes.csv"))
readr::write_csv(as_tibble(activate(g_tbl, edges)), file.path(out_dir, "subgraph_edges.csv"))

message(">>> Done.")

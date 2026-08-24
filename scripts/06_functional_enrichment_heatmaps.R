#!/usr/bin/env Rscript
###############################################################################
# Heatmap & Dendrogram Visualization (looped over enrichment types)
# - Automatic sizing (no manual cell size)
# - Per-cell borders + frames around body/row/column name strips
# - Extra space for axis labels; legends don’t overlap
#
# Inputs (per type): outputs/Enrichment/_summaries/merged_<TYPE>_significant.csv
# Outputs per type  : outputs/Comparisons/<TYPE>/
###############################################################################

source(file.path("scripts", "_config.R"))
assert_packages(c("tidyverse", "ComplexHeatmap", "circlize"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ========================== Helpers ==========================================
make_matrix <- function(data, row_id="Description", col_id="Experiment", value_col="neglog10_padj") {
  if (nrow(data) == 0) return(matrix(0, nrow=0, ncol=0))
  wide <- data %>%
    dplyr::select(all_of(c(row_id, col_id, value_col))) %>%
    group_by(across(all_of(c(row_id, col_id)))) %>%
    summarise(value = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(col_id), values_from = value, values_fill = 0)
  mat <- as.matrix(wide[,-1, drop = FALSE])
  rownames(mat) <- wide[[1]]
  mat
}

col_fun <- colorRamp2(c(0, 5, 10, 20, 50),
                      c("grey90", "lightpink", "red", "darkred", "black"))

# thin border for each cell
.cell_border_fun <- function(j, i, x, y, w, h, fill) {
  grid.rect(x=x, y=y, width=w, height=h, gp=gpar(col="grey70", fill=NA, lwd=0.4))
}

# Space for axis labels (computed from actual text widths)
axis_label_space <- function(mat, row_fs=9, col_fs=10, col_rot=-45) {
  list(
    row_names_gp            = gpar(fontsize=row_fs),
    row_names_max_width     = max_text_width(rownames(mat), gp=gpar(fontsize=row_fs)) + unit(20, "mm"),
    column_names_gp         = gpar(fontsize=col_fs),
    column_names_rot        = col_rot,
    column_names_max_height = max_text_width(colnames(mat), gp=gpar(fontsize=col_fs)) + unit(8, "mm")
  )
}

make_top_anno <- function(mat, exp_info, type_colors,
                          legend_title_fs = 18,
                          legend_label_fs = 14,
                          legend_grid_mm  = 7) {
  exp_names <- colnames(mat)
  exp_types <- exp_info$Type[match(exp_names, exp_info$Experiment)]
  HeatmapAnnotation(
    Type = exp_types,
    col  = list(Type = type_colors),
    annotation_name_side = "right",
    annotation_legend_param = list(
      title_gp    = gpar(fontsize = legend_title_fs, fontface = "bold"),
      labels_gp   = gpar(fontsize = legend_label_fs),
      grid_width  = unit(legend_grid_mm, "mm"),
      grid_height = unit(legend_grid_mm, "mm"),
      direction   = "vertical",     # stack items top-to-bottom
      border      = "black"         # legend box border (optional)
    )
  )
}



# Build Heatmap (no explicit width/height => ComplexHeatmap adapts cells to device)
make_heatmap_with_borders <- function(mat, hm_name,
                                      show_rows = FALSE, row_fs = 8, col_fs = 11,
                                      col_map = col_fun,        # <— was col_fun
                                      top_ann = NULL, col_title = NULL, row_title = NULL,
                                      legend_title_fs = 18,
                                      legend_label_fs = 14,
                                      legend_grid_mm  = 10,
                                      ...) {
  labs <- axis_label_space(mat, row_fs = row_fs, col_fs = col_fs, col_rot = -45)
  ComplexHeatmap::Heatmap(
    mat, name = hm_name, col = col_map,  # <— use col_map here
    cluster_rows = TRUE, cluster_columns = TRUE,
    use_raster = FALSE, layer_fun = .cell_border_fun, border = TRUE,
    show_row_names         = show_rows,
    row_names_gp           = labs$row_names_gp,
    row_names_max_width    = labs$row_names_max_width,
    row_names_side         = "right",
    row_names_centered     = FALSE,
    column_names_gp         = labs$column_names_gp,
    column_names_rot        = labs$column_names_rot,
    column_names_max_height = labs$column_names_max_height,
    heatmap_legend_param = list(
      title_gp    = gpar(fontsize = legend_title_fs, fontface = "bold"),
      labels_gp   = gpar(fontsize = legend_label_fs),
      grid_width  = unit(legend_grid_mm, "mm"),
      grid_height = unit(legend_grid_mm, "mm"),
      direction   = "vertical",
      border      = "black"
    ),
    top_annotation = top_ann,
    column_title   = col_title,
    row_title      = row_title,
    ...
  )
}



# Put ALL, then UP, then DOWN, then Other
get_col_split <- function(mat) {
  ty <- exp_info$Type[match(colnames(mat), exp_info$Experiment)]
  ty[is.na(ty)] <- "Other"
  factor(ty, levels = c("ALL","UP","DOWN","Other"))
}

# Auto page size from matrix aspect ratio (so cells are ~square) unless you pass width/height
# width/height are in inches; if NULL they are computed.
export_and_store <- function(ht, hm_name, filename, outdir,
                             width = NULL, height = NULL, res = 500,
                             mat = NULL,                    # <- NEW (optional)
                             min_width = 10, max_width = 60,
                             base_height = 16,             # used if height is NULL
                             padding_mm = c(12, 28, 20, 20)) {
  
  # ------- derive sensible width/height if not provided -------
  if (is.null(height)) height <- base_height
  if (is.null(width)) {
    if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      aspect <- ncol(mat) / max(1, nrow(mat))   # body aspect ~ columns/rows
      width  <- height * aspect + 6             # +6in for legends/row labels
      width  <- max(min_width, min(width, max_width))
    } else {
      width <- max(min_width, min(base_height, max_width))
    }
  }
  
  draw_once <- function() {
    draw(ht,
         heatmap_legend_side    = "right",
         annotation_legend_side = "right",
         merge_legends          = TRUE,
         padding                = unit(padding_mm, "mm"))
    # optional frames
    try(decorate_heatmap_body(hm_name, { grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
    try(decorate_row_names(hm_name,    { grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
    try(decorate_column_names(hm_name, { grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
  }
  
  # PDF
  pdf(file.path(outdir, paste0(filename, ".pdf")), width = width, height = height)
  draw_once(); dev.off()
  
  # PNG (prefer ragg)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file.path(outdir, paste0(filename, ".png")),
                  width = width * res, height = height * res, units = "px", res = res, background = "white")
    draw_once(); dev.off()
  } else {
    if (!capabilities("cairo")) options(bitmapType = "cairo")
    png(file.path(outdir, paste0(filename, ".png")),
        width = width * res, height = height * res, res = res, type = "cairo")
    draw_once(); dev.off()
  }
  
  invisible(TRUE)
}


# Nicer row titles per enrichment
row_titles <- c(GO_BP="GO BP terms", GO_MF="GO MF terms", GO_CC="GO CC terms",
                KEGG="KEGG pathways", Reactome="Reactome pathways")

# ========================== Loop over enrichment types ========================
ENRICHMENTS <- c("GO_BP","GO_MF","GO_CC","KEGG","Reactome")

for (etype in ENRICHMENTS) {
  infile <- file.path(ENR_DIR, "_summaries", sprintf("merged_%s_significant.csv", etype))
  if (!file.exists(infile)) { message("Skip (missing): ", infile); next }
  
  df <- read.csv(infile, check.names = FALSE)
  if (!all(c("Description","Experiment","p.adjust") %in% names(df))) {
    message("Skip (columns missing) in ", infile); next
  }
  df$neglog10_padj <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  
  # experiment types for annotation
  exp_info <- data.frame(Experiment = unique(df$Experiment)) %>%
    mutate(Type = case_when(
      str_detect(Experiment, "_UP")   ~ "UP",
      str_detect(Experiment, "_DOWN") ~ "DOWN",
      str_detect(Experiment, "_ALL")  ~ "ALL",
      TRUE ~ "Other"
    ))
  type_colors <- c("UP"="firebrick","DOWN"="royalblue","ALL"="darkgreen","Other"="grey40")
  
  outdir <- file.path(CMP_DIR, etype)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Global
  mat_all <- make_matrix(df, row_id="Description", col_id="Experiment")
  top_ann <- make_top_anno(mat_all, exp_info, type_colors)
  hm_name <- "–log10(padj)" #etype
  ht_global <- make_heatmap_with_borders(
    mat_all, hm_name,
    top_ann=top_ann, col_title="Experiments",
    show_rows = TRUE,
    row_fs    = 8,        # keep or adjust
    row_title = row_titles[[etype]] %||% etype,
    column_split = get_col_split(mat_all)
  )
  export_and_store(ht_global, hm_name, sprintf("%s_heatmap_global", etype), outdir,
                   width=40, height=30, mat=mat_all)  # auto-size
  
  # ---- Separate UP / DOWN / ALL
  split_df <- list(
    UP   = df %>% filter(str_detect(Experiment, "_UP")),
    DOWN = df %>% filter(str_detect(Experiment, "_DOWN")),
    ALL  = df %>% filter(str_detect(Experiment, "_ALL"))
  )
  for (part in names(split_df)) {
    dfi <- split_df[[part]]
    if (nrow(dfi) == 0) next
    mati <- make_matrix(dfi)
    topi <- make_top_anno(mati, exp_info, type_colors)
    hm_i <- "–log10(padj)" #paste0(etype, "_", part)
    ht_i <- make_heatmap_with_borders(
      mati, hm_i,top_ann=topi,
      show_rows = TRUE,
      row_fs    = 8,        # keep or adjust
      col_title=paste(part, "genes"),
      column_split = get_col_split(mati)
    )
    export_and_store(ht_i, hm_i, sprintf("%s_heatmap_%s", etype, part), outdir,
                     width=40, height=30, mat=mati)  # auto-size
  }
  
  # ---- Top 20 terms per experiment
  topN <- 10
  df_top <- df %>%
    group_by(Experiment) %>%
    slice_max(order_by = neglog10_padj, n = topN, with_ties = FALSE) %>%
    ungroup()
  if (nrow(df_top)) {
    mat_top <- make_matrix(df_top)
    top_ann2 <- make_top_anno(mat_top, exp_info, type_colors)
    hm_top  <- "–log10(padj)"
    ht_top <- make_heatmap_with_borders(
      mat_top, hm_top, show_rows=TRUE, row_fs=8,
      top_ann=top_ann2,
      col_title = paste("Top", topN, if (etype %in% c("KEGG","Reactome")) "pathways" else "terms", "per Experiment"),
      row_title = row_titles[[etype]] %||% etype,
      column_split = get_col_split(mat_top)
    )
    export_and_store(ht_top, hm_top, sprintf("%s_heatmap_Top20", etype), outdir,
                     width=40, height=30, mat=mat_top)  # auto-size
  }
}

cat("✓ All enrichment heatmaps written under outputs/Comparisons/<TYPE>/\n")

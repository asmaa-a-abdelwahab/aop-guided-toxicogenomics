#!/usr/bin/env Rscript
###############################################################################
# Automated Heatmap generation and export (merged enrichment data)
# Heatmaps + Dendrograms + Combined Multi-panel Report
###############################################################################

source(file.path("scripts", "_config.R"))
assert_packages(c("tidyverse", "ComplexHeatmap", "circlize", "RColorBrewer", "readr"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# ========================== Load & Merge Data ================================
fingerprint_files <- c(
  file.path(AOP_DIR, "AOP_fingerprint_GSE_enriched.csv"),
  file.path(AOP_DIR, "AOP_fingerprint_Literature_enriched.csv")
)
missing_fingerprints <- fingerprint_files[!file.exists(fingerprint_files)]
if (length(missing_fingerprints)) {
  stop(
    "AOP fingerprint tables not found:\n- ",
    paste(missing_fingerprints, collapse = "\n- "),
    "\nRun step 04 to completion before step 05:\n",
    "  Rscript run_pipeline.R --steps=4,5,6",
    call. = FALSE
  )
}
df_gse <- read.csv(fingerprint_files[[1]])
df_lit <- read.csv(fingerprint_files[[2]])
df <- bind_rows(df_gse, df_lit)

# ========================== Output directory =================================
outdir <- file.path(CMP_DIR, "AOP")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ComplexHeatmap text measurement can open R's default file device in a
# non-interactive session. Keep a null device underneath explicit export
# devices so the pipeline never leaves an unintended Rplots.pdf in the root.
grDevices::pdf(NULL)
null_device <- grDevices::dev.cur()

# ========================== Size controls (edit here) ========================
# Per-cell size (mm). Set different values for rectangular cells.
CELL_W_MM      <- 4     # cell width  in mm
CELL_H_MM      <- 24     # cell height in mm

# Optionally constrain total heatmap body size (mm). Set to NA to let cell size rule.
TOTAL_WIDTH_MM  <- 600   # e.g., 260
TOTAL_HEIGHT_MM <- 410   # e.g., 180

# Extra space for names (prevents overlaps). Increase if row/col names are long.
ROW_NAMES_MAX_WIDTH_MM    <- 360
COLUMN_NAMES_MAX_HEIGHT_MM <- 130
TOP_ANNOT_HEIGHT_MM        <- 8   # height of the “Type” annotation bar

# ========================== Helpers =========================================
make_matrix <- function(data, row_id, col_id, value_col = "neglog10_padj") {
  wide <- data %>%
    dplyr::select(all_of(c(row_id, col_id, value_col))) %>%
    group_by(across(all_of(c(row_id, col_id)))) %>%
    summarise(value = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(col_id), values_from = value, values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide[[1]]
  mat
}


col_fun <- colorRamp2(c(0, 5, 10, 20, 50),
                      c("grey90", "lightpink", "red", "darkred", "black"))

# Experiment type annotation
exp_info <- data.frame(Experiment = unique(df$Experiment)) %>%
  mutate(Type = case_when(
    str_detect(Experiment, "_UP")   ~ "UP",
    str_detect(Experiment, "_DOWN") ~ "DOWN",
    str_detect(Experiment, "_ALL")  ~ "ALL",
    TRUE ~ "Other"
  ))

type_colors <- c("UP"="firebrick", "DOWN"="royalblue",
                 "ALL"="darkgreen", "Other"="grey40")

col_anno_fun <- function(mat) {
  exp_names <- colnames(mat)
  exp_types <- exp_info$Type[match(exp_names, exp_info$Experiment)]
  HeatmapAnnotation(
    Type = exp_types,
    col = list(Type = type_colors),
    annotation_name_side = "right",
    height = unit(TOP_ANNOT_HEIGHT_MM, "mm")
  )
}

# ---- safe unit helper -------------------------------------------------------
as_mm_unit <- function(x, default_mm) {
  if (is.null(x)) return(grid::unit(default_mm, "mm"))
  if (inherits(x, "unit")) return(x)
  grid::unit(as.numeric(x), "mm")
}

# ---- body sizing ------------------------------------------------------------
heatmap_dims <- function(mat, cell_w_mm = CELL_W_MM, cell_h_mm = CELL_H_MM,
                         width_mm = TOTAL_WIDTH_MM, height_mm = TOTAL_HEIGHT_MM) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  nr <- nrow(mat); nc <- ncol(mat)
  if (!length(nr) || !length(nc) || nr == 0 || nc == 0)
    return(list(width = grid::unit(0, "mm"), height = grid::unit(0, "mm")))
  nat_w <- nc * as.numeric(cell_w_mm)
  nat_h <- nr * as.numeric(cell_h_mm)
  to_num <- function(x) if (is.null(x) || is.na(x)) NULL else if (inherits(x,"unit")) grid::convertUnit(x,"mm",valueOnly=TRUE) else as.numeric(x)
  W <- to_num(width_mm); H <- to_num(height_mm)
  if (!is.null(W) && !is.null(H)) { w <- W; h <- H
  } else if (!is.null(W))        { s <- W/nat_w; w <- W; h <- nat_h*s
  } else if (!is.null(H))        { s <- H/nat_h; h <- H; w <- nat_w*s
  } else                         { w <- nat_w; h <- nat_h }
  list(width = grid::unit(w, "mm"), height = grid::unit(h, "mm"))
}

# ---- per-cell border --------------------------------------------------------
.cell_border_fun <- function(j, i, x, y, w, h, fill) {
  grid::grid.rect(x=x, y=y, width=w, height=h,
                  gp=grid::gpar(col="grey70", fill=NA, lwd=0.4))
}

# ---- robust builder (accepts ... and filters duplicates) --------------------
make_heatmap_with_borders <- function(
    mat, hm_name,
    col,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = FALSE,
    top_annotation = NULL,
    row_title = NULL, column_title = NULL,
    # label sizing
    row_fs = 9, col_fs = 10, col_rot = 90,
    row_names_max_width = NULL,
    column_names_max_height = NULL,
    ...
) {
  # Defaults based on actual text
  rn_default <- if (!is.null(rownames(mat)) && length(rownames(mat))) {
    ComplexHeatmap::max_text_width(rownames(mat), gp=grid::gpar(fontsize=row_fs)) + grid::unit(6, "mm")
  } else grid::unit(20, "mm")
  cn_default <- if (!is.null(colnames(mat)) && length(colnames(mat))) {
    ComplexHeatmap::max_text_width(colnames(mat), gp=grid::gpar(fontsize=col_fs)) + grid::unit(6, "mm")
  } else grid::unit(12, "mm")
  
  rnw <- if (is.null(row_names_max_width)) rn_default else as_mm_unit(row_names_max_width, 20)
  cnh <- if (is.null(column_names_max_height)) cn_default else as_mm_unit(column_names_max_height, 12)
  
  dims <- heatmap_dims(mat)
  
  base_args <- list(
    mat = mat,
    name = hm_name,
    col  = col,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    use_raster = FALSE,
    layer_fun  = .cell_border_fun,
    border = TRUE,
    width  = dims$width,
    height = dims$height,
    
    show_row_names         = show_row_names,
    row_names_gp           = grid::gpar(fontsize = row_fs),
    row_names_max_width    = rnw,
    column_names_gp        = grid::gpar(fontsize = col_fs),
    column_names_rot       = col_rot,
    column_names_max_height= cnh,
    
    top_annotation = top_annotation,
    row_title = row_title,
    column_title = column_title
  )
  
  dots <- list(...)
  # Remove any duplicates that would clash with our base args
  if (length(dots)) {
    dup_keys <- intersect(names(dots), names(base_args))
    if (length(dup_keys)) dots[dup_keys] <- NULL
  }
  
  do.call(ComplexHeatmap::Heatmap, c(base_args, dots))
}

# Put ALL, then UP, then DOWN, then Other
get_col_split <- function(mat) {
  ty <- exp_info$Type[match(colnames(mat), exp_info$Experiment)]
  ty[is.na(ty)] <- "Other"
  factor(ty, levels = c("ALL","UP","DOWN","Other"))
}


# robust exporter (PDF + PNG w/ downscaling) and frame decorations
export_and_store <- function(ht, hm_name, filename,
                             width = 34, height = 18, res = 500,
                             max_pixels = 80e6,
                             save_svg = TRUE, svg_bg = "white") {
  
  draw_and_box <- function() {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side    = "right",
      annotation_legend_side = "right",
      merge_legends          = TRUE,
      padding                = grid::unit(c(6, 10, 6, 6), "mm")
    )
    try(ComplexHeatmap::decorate_heatmap_body(hm_name, { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 1)) }), silent = TRUE)
    try(ComplexHeatmap::decorate_column_names(hm_name, { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 1)) }), silent = TRUE)
    try(ComplexHeatmap::decorate_annotation("Type",    { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 1)) }), silent = TRUE)
    for (s in 1:50) {
      try(ComplexHeatmap::decorate_row_names(hm_name, slice = s, {
        grid::rect(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   width = unit(1, "npc"), height = unit(1, "npc"),
                   gp = gpar(fill = NA, col = "black", lwd = 0.8, linejoin = "round"))
      }), silent = TRUE)
    }
  }
  
  # ---- PDF ----
  pdf(file.path(outdir, paste0(filename, ".pdf")), width = width, height = height)
  draw_and_box(); dev.off()
  
  # # ---- SVG (vector) ----
  # if (isTRUE(save_svg)) {
  #   svg_path <- file.path(outdir, paste0(filename, ".svg"))
  #   if (requireNamespace("svglite", quietly = TRUE)) {
  #     svglite::svglite(svg_path, width = width, height = height, bg = svg_bg)
  #   } else {
  #     # base device fallback
  #     grDevices::svg(svg_path, width = width, height = height, bg = svg_bg)
  #   }
  #   draw_and_box(); dev.off()
  # }
  
  # ---- PNG (raster, safe downscale) ----
  px_w <- width * res; px_h <- height * res
  scale <- sqrt((px_w * px_h) / max_pixels)
  if (is.finite(scale) && scale > 1) {
    res  <- max(72, round(res/scale))
    px_w <- round(px_w/scale); px_h <- round(px_h/scale)
  }
  png_path <- file.path(outdir, paste0(filename, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_path, width = px_w, height = px_h, units = "px", res = res, background = "white")
    draw_and_box(); dev.off()
  } else {
    if (!capabilities("cairo")) options(bitmapType = "cairo")
    png(png_path, width = px_w, height = px_h, res = res, type = "cairo")
    draw_and_box(); dev.off()
  }
  invisible(ht)
}


plots_list <- list()


topN <- 20
batchN <- 2

# --- helpers: coerce width/height to numeric inches safely -------------------
.to_inches <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (inherits(x, "unit")) return(grid::convertUnit(x, "in", valueOnly = TRUE))
  as.numeric(x)
}

export_and_store <- function(ht, hm_name, filename, outdir,
                             width = 40, height = 30, res = 500,
                             mat = NULL,                    # optional, for autosize logic
                             min_width = 10, max_width = 60,
                             base_height = 16,
                             padding_mm = c(12, 28, 20, 20)) {
  
  # ----- autosize (optional) -------------------------------------------------
  H <- .to_inches(height); W <- .to_inches(width)
  if (is.na(H)) H <- base_height
  if (is.na(W)) {
    if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      aspect <- ncol(mat) / max(1, nrow(mat))   # body aspect ~ columns/rows
      W <- H * aspect + 6                       # + space for legends/labels
      W <- max(min_width, min(W, max_width))
    } else {
      W <- max(min_width, min(base_height, max_width))
    }
  }
  
  # ensure plain numerics for devices
  RES <- as.numeric(res)
  if (!is.finite(W) || !is.finite(H) || !is.finite(RES)) stop("Bad device size/res.")
  
  draw_once <- function() {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side    = "right",
      annotation_legend_side = "right",
      merge_legends          = TRUE,
      padding                = grid::unit(padding_mm, "mm")
    )
    # optional frames
    try(ComplexHeatmap::decorate_heatmap_body(hm_name, { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
    try(ComplexHeatmap::decorate_row_names(hm_name,    { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
    try(ComplexHeatmap::decorate_column_names(hm_name, { grid::rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) }), silent = TRUE)
  }
  
  # PDF
  pdf(file.path(outdir, paste0(filename, ".pdf")), width = W, height = H)
  draw_once(); dev.off()
  
  # PNG (prefer ragg)
  px_w <- as.integer(round(W * RES))
  px_h <- as.integer(round(H * RES))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file.path(outdir, paste0(filename, ".png")),
                  width = px_w, height = px_h, units = "px", res = RES, background = "white")
    draw_once(); dev.off()
  } else {
    if (!capabilities("cairo")) options(bitmapType = "cairo")
    png(file.path(outdir, paste0(filename, ".png")),
        width = px_w, height = px_h, res = RES, type = "cairo")
    draw_once(); dev.off()
  }
  
  # SVG (optional)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(outdir, paste0(filename, ".svg")), width = W, height = H, bg = "white")
    draw_once(); dev.off()
  }

  plots_list[[filename]] <<- ht
  
  invisible(TRUE)
}



# =================== Top-N in rolling chunks (KEs & AOPs) ====================

# Helper: rank per experiment and split into N-sized batches
rank_in_chunks <- function(df, id_col, value_col = "neglog10_padj", topN = 20) {
  stopifnot(all(c("Experiment", id_col, value_col) %in% names(df)))
  df %>%
    group_by(Experiment, .data[[id_col]]) %>%
    summarise(!!value_col := max(.data[[value_col]], na.rm = TRUE), .groups = "drop_last") %>%
    arrange(desc(.data[[value_col]]), .by_group = TRUE) %>%
    mutate(rank  = row_number(),
           batch = ceiling(rank / topN)) %>%
    ungroup()
}

# Build one deterministic annotation per KE/AOP. Missing annotations fall back
# to "Uncategorized"; ties are resolved alphabetically after frequency.
dominant_group <- function(data, id_col, group_col) {
  stopifnot(id_col %in% names(data))

  ids <- tibble::tibble(id = as.character(data[[id_col]])) %>%
    filter(!is.na(id), nzchar(id)) %>%
    distinct()

  if (!group_col %in% names(data)) {
    return(mutate(ids, Group = "Uncategorized"))
  }

  counts <- tibble::tibble(
    id = as.character(data[[id_col]]),
    Group = as.character(data[[group_col]])
  ) %>%
    filter(!is.na(id), nzchar(id), !is.na(Group), nzchar(Group)) %>%
    count(id, Group, name = ".n") %>%
    arrange(id, desc(.n), Group) %>%
    group_by(id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(id, Group)

  ids %>%
    left_join(counts, by = "id", relationship = "one-to-one") %>%
    mutate(Group = coalesce(Group, "Uncategorized"))
}

# Render all batches for a given ranked table
render_top_batches <- function(ranked_tbl, id_col, prefix, outdir,
                               title_entity, topN = 20, batchN = 5) {
  batches <- sort(unique(ranked_tbl$batch))
  batches <- batches[seq_len(min(batchN, length(batches)))]
  if (length(batches) == 0) return(invisible(NULL))
  
  for (b in batches) {
    df_chunk <- ranked_tbl %>% filter(batch == b) %>% ungroup()
    if (nrow(df_chunk) == 0) next
    
    # Make matrix from union of per-experiment TopN of this batch
    mat <- make_matrix(df_chunk, row_id = id_col, col_id = "Experiment", value_col = "neglog10_padj")
    if (nrow(mat) == 0 || ncol(mat) == 0) next
    
    # Annotation + heatmap
    top_ann <- col_anno_fun(mat)
    hm_name <- "–log10(padj)"
    ht <- make_heatmap_with_borders(
      mat, hm_name,
      col = col_fun,                     # << use your col_fun
      cluster_rows = TRUE, cluster_columns = TRUE,
      show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
      top_annotation = top_ann,
      column_title = paste0(
        if (b == 1) "Top " else paste0(b, "nd Top "),  # quick label; tweak if you want 2nd/3rd suffixes perfect
        topN, " ", title_entity, " per Experiment (Batch ", b, ")"
      ),
      row_title = title_entity,
      column_split = get_col_split(mat)
    )
    
    # File names + CSV
    file_stub <- sprintf("%s_Top%02d_Batch%02d", prefix, topN, b)
    readr::write_csv(
      df_chunk %>% dplyr::select(Experiment, !!rlang::sym(id_col), neglog10_padj, rank, batch),
      file.path(outdir, paste0(file_stub, ".csv"))
    )
    export_and_store(ht, hm_name, file_stub, outdir, width = 40, height = 30, mat = mat)
  }
  invisible(TRUE)
}

# ------------------------------ KEs ------------------------------------------

# Rank KEs in rolling Top-N chunks
ranked_ke <- df %>%
  filter(!is.na(Ke_description)) %>%
  rank_in_chunks(id_col = "Ke_description", value_col = "neglog10_padj", topN = topN)

# (Optional) also keep the single "Top N" products you already had:
df_top_ke <- ranked_ke %>% filter(batch == 1) %>% dplyr::select(-rank, -batch)
readr::write_csv(df_top_ke, file.path(outdir, sprintf("Top%02d_KEs_per_experiment.csv", topN)))

# Render *all* batches for KEs
render_top_batches(
  ranked_tbl   = ranked_ke,
  id_col       = "Ke_description",
  prefix       = "KEs",
  outdir       = outdir,
  title_entity = "Key Events",
  topN         = topN,
  batchN       = batchN
)

# ------------------------------ AOPs -----------------------------------------
ranked_aop <- df %>%
  filter(!is.na(AOP_name)) %>%
  rank_in_chunks(id_col = "AOP_name", value_col = "neglog10_padj", topN = topN)

df_top_aop <- ranked_aop %>% filter(batch == 1) %>% dplyr::select(-rank, -batch)
readr::write_csv(df_top_aop, file.path(outdir, sprintf("Top%02d_AOPs_per_experiment.csv", topN)))

render_top_batches(
  ranked_tbl   = ranked_aop,
  id_col       = "AOP_name",
  prefix       = "AOPs",
  outdir       = outdir,
  title_entity = "AOPs",
  topN         = topN,
  batchN       = batchN
)



# ---------- helper to render grouped Top-N batches (KEs or AOPs) ------------
render_grouped_top_batches <- function(ranked_tbl, id_col, group_map, group_label,
                                       prefix, outdir, topN = 20, batchN = 5) {
  # ensure one-to-one map
  group_map <- group_map %>% distinct(id, .keep_all = TRUE)
  
  batches <- sort(unique(ranked_tbl$batch))
  batches <- batches[seq_len(min(batchN, length(batches)))]
  
  for (b in batches) {
    df_chunk <- ranked_tbl %>% filter(batch == b)
    if (nrow(df_chunk) == 0) next
    
    # attach group (SSbD/Organ/Endpoint) to the item id
    df_chunk2 <- df_chunk %>%
      left_join(group_map, by = setNames("id", id_col), relationship = "many-to-one") %>%
      mutate(Group = dplyr::coalesce(Group, "Uncategorized")) %>%
      distinct(Experiment, .data[[id_col]], neglog10_padj, rank, batch, Group)
    
    # build matrix + row splits
    mat <- make_matrix(df_chunk2, row_id = id_col, col_id = "Experiment",
                       value_col = "neglog10_padj")
    if (nrow(mat) == 0 || ncol(mat) == 0) next
    
    # row split factor in the same order across batches
    rg_tbl <- df_chunk2 %>% distinct(!!rlang::sym(id_col), Group)
    row_groups <- rg_tbl$Group[match(rownames(mat), rg_tbl[[id_col]])]
    row_groups[is.na(row_groups) | row_groups == ""] <- "Uncategorized"
    row_groups <- factor(row_groups, levels = sort(unique(row_groups)))
    
    # right-side labeled blocks (labels horizontal)
    right_blocks <- rowAnnotation(
      Group = anno_block(
        gp         = gpar(fill = NA, col = "black", lwd = 1.2),
        labels     = levels(row_groups),
        labels_gp  = gpar(fontsize = 9, fontface = "bold"),
        labels_rot = 0
      ),
      width = unit(70, "mm")
    )
    
    # column annotation + heatmap
    top_ann <- col_anno_fun(mat)
    hm_name <- "–log10(padj)"
    
    # a bit of space for row names inside their box
    row_name_space <- max_text_width(rownames(mat), gp = gpar(fontsize = 9)) + unit(6, "mm")
    
    ht <- make_heatmap_with_borders(
      mat, hm_name,
      col       = col_fun,                         # use your palette
      show_row_names = TRUE, row_names_gp = gpar(fontsize = 9),
      row_names_max_width = row_name_space,
      row_split  = row_groups,
      right_annotation = right_blocks,             # groups OUTSIDE the box
      top_annotation = top_ann,
      column_title = sprintf("Batch %d (Top-%d window) — grouped by %s", b, topN, group_label),
      row_title    = group_label,
      column_split = get_col_split(mat)            # ALL / UP / DOWN / Other together
    )
    
    # write CSV for traceability
    stub <- sprintf("%s_GroupedBy_%s_Top%02d_Batch%02d",
                    prefix, gsub("[^A-Za-z0-9]+","_", tolower(group_label)), topN, b)
    readr::write_csv(df_chunk2, file.path(outdir, paste0(stub, ".csv")))
    
    # export (PDF/PNG/SVG if your export_and_store writes SVG too)
    export_and_store(ht, hm_name, stub, outdir, width = 40, height = 28, mat = mat)
  }
  invisible(TRUE)
}

# -------------------- KEs: grouped batches by SSbD / Organ / Endpoint -------
ke_maps <- list(
  SSbD = dominant_group(df, "Ke_description", "SSbD_category"),
  Organ = dominant_group(df, "Ke_description", "Organ"),
  Endpoint = dominant_group(df, "Ke_description", "Endpoint")
)

render_grouped_top_batches(ranked_ke, "Ke_description", ke_maps$SSbD,     "SSbD category",
                           prefix = "KEs", outdir = outdir, topN = topN, batchN = batchN)
render_grouped_top_batches(ranked_ke, "Ke_description", ke_maps$Organ,    "Organ",
                           prefix = "KEs", outdir = outdir, topN = topN, batchN = batchN)
render_grouped_top_batches(ranked_ke, "Ke_description", ke_maps$Endpoint, "Endpoint",
                           prefix = "KEs", outdir = outdir, topN = topN, batchN = batchN)

# -------------------- AOPs: grouped batches by SSbD / Organ / Endpoint ------
aop_maps <- list(
  SSbD = dominant_group(df, "AOP_name", "SSbD_category"),
  Organ = dominant_group(df, "AOP_name", "Organ"),
  Endpoint = dominant_group(df, "AOP_name", "Endpoint")
)

render_grouped_top_batches(ranked_aop, "AOP_name", aop_maps$SSbD,     "SSbD category",
                           prefix = "AOPs", outdir = outdir, topN = topN, batchN = batchN)
render_grouped_top_batches(ranked_aop, "AOP_name", aop_maps$Organ,    "Organ",
                           prefix = "AOPs", outdir = outdir, topN = topN, batchN = batchN)
render_grouped_top_batches(ranked_aop, "AOP_name", aop_maps$Endpoint, "Endpoint",
                           prefix = "AOPs", outdir = outdir, topN = topN, batchN = batchN)


###############################################################################
# Combined multi-panel PDF report
###############################################################################
pdf(file.path(outdir, "heatmap_report_combined.pdf"), width = 28, height = 16)
for (nm in names(plots_list)) {
  grid.newpage()
  grid.text(nm, x = 0.5, y = 0.96, gp = gpar(fontsize = 16, fontface = "bold"))
  draw(plots_list[[nm]], heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(6, 10, 6, 6), "mm"))
}
invisible(dev.off())

if (identical(grDevices::dev.cur(), null_device)) {
  invisible(grDevices::dev.off())
}

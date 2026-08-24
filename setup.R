#!/usr/bin/env Rscript

# Bootstrap a reproducible renv environment for the pipeline.
project <- normalizePath(".", winslash = "/", mustWork = TRUE)
cran_repo <- "https://cloud.r-project.org"

if (!file.exists(file.path(project, "scripts", "_config.R"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

options(
  repos = c(CRAN = cran_repo),
  Ncpus = max(1L, parallel::detectCores() - 1L),
  timeout = max(600, getOption("timeout", 60))
)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = cran_repo)
}

renv::consent(provided = TRUE)
if (!file.exists(file.path(project, "renv", "activate.R"))) {
  renv::init(project = project, bare = TRUE, restart = FALSE)
} else {
  renv::load(project = project)
}

# R 4.4 is paired with Bioconductor 3.20 for this project.
if (getRversion() < "4.4" || getRversion() >= "4.5") {
  warning(
    "This project is tested with R 4.4.x. Current R version: ",
    as.character(getRversion())
  )
}
renv::settings$bioconductor.version("3.20", persist = TRUE)

cran_packages <- c(
  "BiocManager", "circlize", "dplyr", "forcats", "fs", "ggraph",
  "ggplot2", "htmltools", "htmlwidgets", "igraph", "jsonlite",
  "matrixStats", "msigdbr", "openxlsx", "optparse", "pheatmap", "purrr",
  "RColorBrewer", "ragg", "readr", "readxl", "scales", "stringr",
  "svglite", "tibble", "tidygraph", "tidyr", "tidyverse", "visNetwork"
)

bioconductor_packages <- c(
  "AnnotationDbi", "BiocParallel", "biomaRt", "clusterProfiler", "ComplexHeatmap",
  "DESeq2", "dorothea", "enrichplot", "fgsea", "org.Hs.eg.db",
  "ReactomePA", "vsn"
)

specs <- c(
  paste0("cran::", cran_packages),
  paste0("bioc::", bioconductor_packages),
  "kassambara/easyGgplot2",
  "fhaive/AOPfingerprintR"
)

message("Installing/restoring pipeline dependencies...")
renv::install(specs)
renv::snapshot(project = project, prompt = FALSE)
renv::status(project = project)

message("Setup complete. Commit renv.lock and the renv bootstrap files to Git.")

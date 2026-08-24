# Shared configuration for the AOP-guided toxicogenomics pipeline.

options(stringsAsFactors = FALSE)
set.seed(123)

ROOT <- normalizePath(
  Sys.getenv("AOP_PIPELINE_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)

if (!file.exists(file.path(ROOT, "scripts", "_config.R"))) {
  stop(
    "Run pipeline scripts from the repository root, or set AOP_PIPELINE_ROOT ",
    "to the repository path.",
    call. = FALSE
  )
}

slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

path_setting <- function(name, default) {
  configured <- trimws(Sys.getenv(name, unset = ""))
  value <- if (nzchar(configured)) configured else default
  if (!grepl("^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)", value)) {
    value <- file.path(ROOT, value)
  }
  normalizePath(value, winslash = "/", mustWork = FALSE)
}

csv_setting <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  result <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  result[nzchar(result)]
}

numeric_setting <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) stop(name, " must be numeric.", call. = FALSE)
  value
}

# Inputs
COUNTS_FILE <- path_setting("AOP_COUNTS_FILE", "data/GSE200036_counts.tsv")
META_FILE   <- path_setting("AOP_METADATA_FILE", "data/GSE200036_metadata.csv")
LIT_FILE    <- path_setting("AOP_LITERATURE_FILE", "data/Literature_DEGs.csv")
BACKGROUND_GSE_FILE <- path_setting("AOP_BACKGROUND_GSE_FILE", "data/background_genes.xlsx")
BACKGROUND_LIT_FILE <- path_setting("AOP_BACKGROUND_LIT_FILE", "data/background_genes_GRCh38.xlsx")

# Outputs
CASE_STUDY_ID <- trimws(Sys.getenv("AOP_CASE_STUDY_ID", unset = ""))
default_output <- if (nzchar(CASE_STUDY_ID)) {
  file.path("outputs", slug(CASE_STUDY_ID))
} else {
  "outputs"
}
OUT_DIR <- path_setting("AOP_OUTPUT_DIR", default_output)
DE_DIR  <- file.path(OUT_DIR, "DEG_tables")
ENR_DIR <- file.path(OUT_DIR, "Enrichment")
CMP_DIR <- file.path(OUT_DIR, "Comparisons")
AOP_DIR <- file.path(OUT_DIR, "AOP")

for (path in c(OUT_DIR, DE_DIR, ENR_DIR, CMP_DIR, AOP_DIR)) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Analysis settings
RUN_ON          <- toupper(csv_setting("AOP_RUN_ON", c("UP", "DOWN", "ALL")))
CHEMS           <- csv_setting("AOP_CHEMICALS", c("Ta4C3", "Mo2Ti2C3", "Nb4C3", "LPS", "ConA"))
CONTROL_ALIASES <- csv_setting("AOP_CONTROL_ALIASES", c("CTRL", "Control", "Medium"))
FDR_CUTOFF      <- numeric_setting("AOP_FDR_CUTOFF", 0.05)
LFC_MIN         <- numeric_setting("AOP_LFC_MIN", 1)
TIMEPOINT       <- numeric_setting("AOP_TIMEPOINT", 24)

if (!length(CHEMS)) stop("AOP_CHEMICALS must contain at least one condition.", call. = FALSE)
if (!length(CONTROL_ALIASES)) stop("AOP_CONTROL_ALIASES must not be empty.", call. = FALSE)
if (!all(RUN_ON %in% c("UP", "DOWN", "ALL"))) {
  stop("AOP_RUN_ON may contain only UP, DOWN, and ALL.", call. = FALSE)
}
if (LFC_MIN < 0) stop("AOP_LFC_MIN must be non-negative.", call. = FALSE)
if (FDR_CUTOFF <= 0 || FDR_CUTOFF > 1) {
  stop("AOP_FDR_CUTOFF must be in (0, 1].", call. = FALSE)
}

strip_v <- function(x) sub("\\..*$", "", x)

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Run `Rscript setup.R` from the repository root.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

assert_files <- function(paths, label = "Required input files") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      label, " not found:\n- ", paste(missing, collapse = "\n- "),
      "\nSee data/README.md for the expected inputs.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

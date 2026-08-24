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

# Inputs
COUNTS_FILE <- file.path(ROOT, "data", "GSE200036_counts.tsv")
META_FILE   <- file.path(ROOT, "data", "GSE200036_metadata.csv")
LIT_FILE    <- file.path(ROOT, "data", "Literature_DEGs.csv")

# Outputs
OUT_DIR <- file.path(ROOT, "outputs")
DE_DIR  <- file.path(OUT_DIR, "DEG_tables")
ENR_DIR <- file.path(OUT_DIR, "Enrichment")
CMP_DIR <- file.path(OUT_DIR, "Comparisons")
AOP_DIR <- file.path(OUT_DIR, "AOP")

for (path in c(OUT_DIR, DE_DIR, ENR_DIR, CMP_DIR, AOP_DIR)) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Analysis settings
RUN_ON     <- c("UP", "DOWN", "ALL")
CHEMS      <- c("Ta4C3", "Mo2Ti2C3", "Nb4C3", "LPS", "ConA")
FDR_CUTOFF <- 0.05
LFC_MIN    <- 1

strip_v <- function(x) sub("\\..*$", "", x)
slug    <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

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

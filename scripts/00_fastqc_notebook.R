#!/usr/bin/env Rscript
# Execute the portable FASTQ QC notebook as an optional upstream stage.

source(file.path("scripts", "_config.R"))

notebook <- path_setting("AOP_FASTQC_NOTEBOOK", "notebooks/Fastq_QC.ipynb")
fastq_dir <- path_setting("AOP_FASTQ_DIR", "data/fastq")
fastqc_out <- path_setting("AOP_FASTQC_OUTPUT_DIR", file.path(OUT_DIR, "FastQC"))
threads <- as.integer(numeric_setting("AOP_FASTQC_THREADS", 4))

assert_files(notebook, "FastQC notebook")
if (!dir.exists(fastq_dir)) {
  stop(
    "FASTQ directory not found: ", fastq_dir,
    "\nSet AOP_FASTQ_DIR to the folder containing *.fastq.gz files.",
    call. = FALSE
  )
}

fastq_files <- unique(unlist(lapply(
  c("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"),
  function(pattern) Sys.glob(file.path(fastq_dir, pattern))
)))
if (!length(fastq_files)) {
  stop(
    "No FASTQ files were found in ", fastq_dir,
    ". This optional stage requires locally available reads.",
    call. = FALSE
  )
}

jupyter <- Sys.which("jupyter")
if (!nzchar(jupyter)) {
  stop(
    "Jupyter was not found on PATH. Install Jupyter/nbconvert or execute ",
    notebook, " manually after setting AOP_FASTQ_DIR.",
    call. = FALSE
  )
}

dir.create(fastqc_out, recursive = TRUE, showWarnings = FALSE)
executed_name <- "Fastq_QC_executed.ipynb"
args <- c(
  "nbconvert", "--to", "notebook", "--execute",
  shQuote(notebook),
  "--output", shQuote(executed_name),
  "--output-dir", shQuote(fastqc_out),
  "--ExecutePreprocessor.timeout=-1"
)
env <- c(
  sprintf("AOP_FASTQ_DIR=%s", fastq_dir),
  sprintf("AOP_FASTQC_OUTPUT_DIR=%s", fastqc_out),
  sprintf("AOP_FASTQC_THREADS=%d", threads)
)

message(sprintf("Executing FastQC notebook for %d FASTQ file(s).", length(fastq_files)))
status <- system2(jupyter, args = args, env = env, stdout = "", stderr = "")
if (!identical(status, 0L)) {
  stop("FastQC notebook execution failed (exit status ", status, ").", call. = FALSE)
}

cat("FastQC notebook complete:", normalizePath(fastqc_out, winslash = "/"), "\n")

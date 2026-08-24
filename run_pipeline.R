#!/usr/bin/env Rscript

# Run the numbered pipeline from the repository root.
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
runner_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "run_pipeline.R"
project_root <- normalizePath(dirname(runner_path), winslash = "/", mustWork = TRUE)
setwd(project_root)
Sys.setenv(AOP_PIPELINE_ROOT = project_root)

args <- commandArgs(trailingOnly = TRUE)

steps <- c(
  `1` = "scripts/01_deseq2_qc.R",
  `2` = "scripts/02_tf_target_expansion.R",
  `3` = "scripts/03_functional_enrichment.R",
  `4` = "scripts/04_aop_fingerprint.R",
  `5` = "scripts/05_aop_heatmaps.R",
  `6` = "scripts/06_functional_enrichment_heatmaps.R",
  `7` = "scripts/07_aop_network.R"
)

usage <- function() {
  cat(
    "Usage: Rscript run_pipeline.R [--steps=1,2,3,4,5,6] [--dry-run]\n\n",
    "Steps 1-6 run by default. Step 7 is optional because it requires the\n",
    "AOP-Wiki TSV exports described in data/README.md.\n",
    sep = ""
  )
}

if (any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

steps_arg <- grep("^--steps=", args, value = TRUE)
selected <- if (length(steps_arg)) {
  strsplit(sub("^--steps=", "", steps_arg[[1]]), ",", fixed = TRUE)[[1]]
} else {
  as.character(1:6)
}
selected <- trimws(selected)

unknown <- setdiff(selected, names(steps))
if (length(unknown)) {
  stop("Unknown step(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

dry_run <- "--dry-run" %in% args
rscript <- file.path(R.home("bin"), "Rscript")

cat("Project root: ", project_root, "\n", sep = "")
cat("Selected steps: ", paste(selected, collapse = ", "), "\n", sep = "")

for (id in selected) {
  script <- normalizePath(steps[[id]], winslash = "/", mustWork = TRUE)
  cat(sprintf("\n[%s/%s] %s\n", match(id, selected), length(selected), basename(script)))
  if (dry_run) next

  status <- system2(rscript, args = shQuote(script), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Pipeline stopped at step ", id, " (exit status ", status, ").", call. = FALSE)
  }
}

if (dry_run) {
  cat("\nDry run complete; no analysis scripts were executed.\n")
} else {
  cat("\nPipeline completed successfully.\n")
}

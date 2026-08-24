# AOP-Guided Toxicogenomics Pipeline

An analysis pipeline that combines transcriptomic differential expression and
literature-derived gene sets with functional enrichment, Key Event enrichment,
Adverse Outcome Pathway (AOP) fingerprinting, heatmaps, and network
visualization.

![AOP-guided toxicogenomics workflow](assets/workflow.png)

## Workflow

0. Validate the count matrix, metadata, sample matching, library sizes, and QC.
1. Run DESeq2 quality control and differential expression.
2. Expand literature-derived gene sets with DoRothEA transcription-factor
   targets at AB and ABC confidence levels.
3. Run GO, KEGG, and Reactome enrichment.
4. Enrich Key Events and create AOP fingerprints.
5. Render cross-experiment AOP heatmaps.
6. Render functional-enrichment heatmaps.
7. Optionally build an AOP-Wiki subnetwork in R or Python.
8. Optionally project literature UP/DOWN signatures onto each transcriptomic
   contrast with preranked FGSEA and render an NES heatmap.
9. Summarize cross-condition DEG burden, overlap, effect-size concordance, and
   recurrent genes with publication-ready figures.
10. Integrate recurrent functional terms, AOP profiles, TF-expansion effects,
    and (when available) literature-signature projection results.

## Repository layout

```text
.
├── assets/                 # Documentation images
├── data/                   # Local inputs (ignored except for documentation)
├── outputs/                # Generated results (ignored)
├── scripts/                # Numbered analysis steps and shared configuration
├── run_pipeline.R          # Sequential pipeline runner
├── setup.R                 # renv dependency bootstrap
└── requirements.txt        # Optional Python network dependencies
```

## Requirements

- R 4.4.x and Bioconductor 3.20.
- System libraries needed by common R graphics/network packages. On
  Debian/Ubuntu these commonly include `libcurl4-openssl-dev`,
  `libfontconfig1-dev`, `libglpk-dev`, and `libssl-dev`.
- Python 3.10 or newer only if using the Python network visualizer.

## Setup

Clone the repository, open a terminal at its root, add the inputs described in
[`data/README.md`](data/README.md), and create the R environment:

```bash
Rscript setup.R
```

The setup script initializes `renv`, installs the CRAN/Bioconductor packages and
`fhaive/AOPfingerprintR`, and writes `renv.lock`. Commit the generated lockfile
and `renv` bootstrap files so collaborators restore the same environment.

For the optional Python implementation of step 7:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

## Run

Run validation plus steps 1-6 in order:

```bash
Rscript run_pipeline.R
```

Select specific steps or preview the execution order:

```bash
Rscript run_pipeline.R --steps=00,1,3,6
Rscript run_pipeline.R --steps=00,1,2,3,4,5,6,7,8,9,10 --dry-run
```

Each numbered script can also be run directly from the repository root. Step 02
and both step-07 scripts expose command-line options:

```bash
Rscript scripts/02_tf_target_expansion.R --help
Rscript scripts/07_aop_network.R --help
python scripts/07_aop_network.py --help
Rscript run_pipeline.R --steps=8,9,10
```

Generated tables, figures, serialized enrichment objects, and interactive HTML
files are written under `outputs/` and are intentionally excluded from Git.
Steps 9 and 10 write PNG figures, editable vector PDFs, plot-source CSV files,
and data-driven Markdown summaries under `outputs/Interpretation/`.

## Using the pipeline for another case study

The numbered stages do not assume a fixed number of exposed conditions. The
condition list, paths, thresholds, timepoint, and output location are supplied
by `scripts/_config.R` and can be overridden with environment variables. Use a
separate `AOP_CASE_STUDY_ID` for each study to prevent results from different
datasets being mixed.

1. Prepare a raw integer count matrix and sample metadata using the schema in
   [`data/README.md`](data/README.md). Metadata must contain `sample_name` and
   `condition`; count columns and metadata sample names must match exactly.
2. Prepare the literature gene table and two analysis-background workbooks.
   Literature chemical names should be represented consistently in the
   `Chemical` column and `AOP_CHEMICALS` setting.
3. Set the case-study parameters below, or edit their defaults in
   `scripts/_config.R`. Relative paths are resolved from the repository root.
4. Run step 00 first and review `Data_validation/QC_report.md`, sample matching,
   replicate counts, library sizes, PCA, and identifier mapping before running
   downstream stages.
5. Run the core workflow with `Rscript run_pipeline.R`. Then run steps 8-10 for
   interpretation. Step 7 is independent and only runs when AOP-Wiki exports
   are available.
6. Review all automatically ranked results against dose, exposure duration,
   biological replicates, and the experimental design before publication.

| Environment variable | Purpose | Default |
|---|---|---|
| `AOP_CASE_STUDY_ID` | Writes to `outputs/<case-study>/` | empty (uses `outputs/`) |
| `AOP_OUTPUT_DIR` | Explicit output directory; overrides the case-study directory | `outputs/` |
| `AOP_COUNTS_FILE` | Raw count matrix | `data/GSE200036_counts.tsv` |
| `AOP_METADATA_FILE` | Sample metadata | `data/GSE200036_metadata.csv` |
| `AOP_LITERATURE_FILE` | Literature gene table | `data/Literature_DEGs.csv` |
| `AOP_BACKGROUND_GSE_FILE` | Transcriptomic gene universe | `data/background_genes.xlsx` |
| `AOP_BACKGROUND_LIT_FILE` | Literature gene universe | `data/background_genes_GRCh38.xlsx` |
| `AOP_CHEMICALS` | Comma-separated exposed conditions | current five-condition study |
| `AOP_CONTROL_ALIASES` | Comma-separated labels normalized to `CTRL` | `CTRL,Control,Medium` |
| `AOP_RUN_ON` | Gene-set directions | `UP,DOWN,ALL` |
| `AOP_FDR_CUTOFF` | Multiple-testing threshold | `0.05` |
| `AOP_LFC_MIN` | Absolute log2-fold-change threshold | `1` |
| `AOP_TIMEPOINT` | Timepoint attached to downstream gene lists | `24` |
| `AOP_PUBLICATION_TOP_N` | Terms shown per functional/annotation panel in step 10 | `10` |
| `AOP_PUBLICATION_TOP_AOPS` | Recurrent AOPs shown in the step-10 landscape | `25` |
| `AOP_PUBLICATION_TOP_GENES` | Recurrent genes shown in the step-09 heatmap | `50` |
| `AOP_PUBLICATION_DPI` | Raster figure resolution for steps 9-10 | `600` |

For example, in a POSIX shell:

```bash
export AOP_CASE_STUDY_ID="study_b"
export AOP_COUNTS_FILE="data/study_b_counts.tsv"
export AOP_METADATA_FILE="data/study_b_metadata.csv"
export AOP_LITERATURE_FILE="data/study_b_literature.csv"
export AOP_CHEMICALS="CompoundA,CompoundB,PositiveControl"
export AOP_TIMEPOINT="48"
Rscript run_pipeline.R --steps=00,1,2,3,4,5,6,8,9,10
```

In PowerShell, set the same values with syntax such as
`$env:AOP_CASE_STUDY_ID = "study_b"`.

Step dependencies are intentionally explicit:

| Step | Requires |
|---|---|
| 00-01 | Counts and metadata |
| 02 | Literature gene table |
| 03 | Outputs from 01 and 02 |
| 04 | Outputs from 01/02 and both background workbooks |
| 05 | Output from 04 |
| 06 | Output from 03 |
| 07 | Output from 04 plus AOP-Wiki TSV exports |
| 08 | Ranked output from 01 plus directional literature genes |
| 09 | At least two all-gene contrast tables from 01 |
| 10 | Outputs from 02-04; it also uses output from 08 when present |

The current implementation is for human bulk RNA-seq: it uses human Ensembl
IDs, `org.Hs.eg.db`, human DoRothEA, KEGG `hsa`, Reactome `human`, and the human
AOP fingerprint resources. A non-human case study requires coordinated changes
to those annotation resources and organism arguments. The DESeq2 model is
`~ condition`; studies with batches, pairing, donor effects, or interactions
must extend the design formula in validation and differential-expression steps
rather than treating those variables as condition labels.

## Reproducibility and publication notes

- The random seed and shared paths are defined in `scripts/_config.R`.
- Raw data and generated outputs are not tracked. Document stable accessions,
  checksums, preprocessing, and download instructions for every released input.
- Review the chemical list and thresholds in `scripts/_config.R` before running
  a different dataset.
- The included GitHub Actions workflow checks R and Python syntax without
  downloading the full scientific software stack.

## License

No software license has been selected. Add a `LICENSE` file before public release
if you want others to be able to reuse or modify the code under stated terms.

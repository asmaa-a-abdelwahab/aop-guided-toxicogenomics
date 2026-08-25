# AOP-Guided Toxicogenomics Pipeline

An analysis pipeline that combines transcriptomic differential expression and
literature-derived gene sets with functional enrichment, Key Event enrichment,
Adverse Outcome Pathway (AOP) fingerprinting, heatmaps, and network
visualization.

![AOP-guided toxicogenomics workflow](assets/workflow.png)

## Workflow

- **FASTQ (optional):** execute the portable FastQC/MultiQC notebook when local
  reads are available. This documents upstream sequence quality; the deposited
  raw-count matrix remains the input to the statistical workflow.

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
11. Test KE/AOP specificity against branch-specific, size-matched random gene
    sets and calculate empirical P values.
12. Assemble manuscript figures as 600-dpi PNG and editable vector PDF files
    from pipeline outputs.

## Repository layout

```text
.
├── assets/                 # Documentation images
├── data/                   # Local inputs (ignored except for documentation)
├── outputs/                # Generated results (ignored)
├── notebooks/              # Portable optional FASTQ QC notebook
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
- FastQC plus Jupyter/nbconvert for the optional `fastqc` stage; MultiQC is
  optional but recommended. These tools are not needed when analysis begins
  from an already quality-controlled count matrix.

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
Rscript run_pipeline.R --steps=fastqc,00,1,2,3,4,5,6,7,8,9,10,11,12 --dry-run
```

Each numbered script can also be run directly from the repository root. Step 02
and both step-07 scripts expose command-line options:

```bash
Rscript scripts/02_tf_target_expansion.R --help
Rscript scripts/07_aop_network.R --help
python scripts/07_aop_network.py --help
Rscript run_pipeline.R --steps=8,9,10,11,12
```

Generated tables, figures, serialized enrichment objects, and interactive HTML
files are written under `outputs/` and are intentionally excluded from Git.
Steps 9 and 10 write PNG figures, editable vector PDFs, plot-source CSV files,
and data-driven Markdown summaries under `outputs/Interpretation/`.
Step 11 writes observed, null, term-level, and provenance tables under
`outputs/Specificity_null/`. Step 12 writes the complete publication figure set
and `publication_figure_manifest.csv` under `outputs/Publication_figures/`.

### Optional FASTQ quality control

The repository includes [`notebooks/Fastq_QC.ipynb`](notebooks/Fastq_QC.ipynb),
which replaces the original machine-specific notebook with an
environment-driven workflow. Put local reads in `data/fastq/` or set
`AOP_FASTQ_DIR`, then run:

```bash
Rscript run_pipeline.R --steps=fastqc
```

The notebook discovers FASTQ/FQ files, runs FastQC, optionally aggregates the
reports with MultiQC, and records an executed notebook, module summary, and run
manifest under `outputs/FastQC/`. FASTQ quality control is deliberately an
optional upstream stage because public read files are large and the core case
study starts from the deposited count matrix.

## Literature evidence selection and provenance

The case-study search was exported from
[PubMed ChemInsight](https://pubmed-cheminsight.edelweissconnect.com/) on 6
September 2025. The unmodified export and the curated evidence workbook are
retained under `manuscript/`. Twenty-seven returned article records were
screened for organism, assay, access to gene-level results, and compatibility
with the human transcriptomic branch. Four publication-based omics candidates
and one public GEO study were documented; the quantitative literature branch
uses the two human studies for which directional gene-level signatures could be
curated (Ti3C2Tx and HN-Ti3C2). The workbook records exclusions and source-study
metadata rather than treating highlighted search results as automatically
eligible evidence.

For another case study, archive the original search export, preserve a screening
decision for every record, and create one row per source in `Study_Metadata`.
Curated gene rows should retain source DOI/accession, material, organism, model,
exposure, direction, identifier type, extraction route, and any transformation.
Do not combine species or assay types merely to increase the apparent evidence
base.

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
5. Run the core workflow with `Rscript run_pipeline.R`. Then run steps 8-12 for
   interpretation, specificity testing, and publication figures. Step 7 is
   independent and only runs when AOP-Wiki exports are available.
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
| `AOP_FASTQ_DIR` | Local FASTQ/FQ directory for the optional notebook stage | `data/fastq/` |
| `AOP_FASTQC_OUTPUT_DIR` | FastQC/MultiQC output directory | `outputs/FastQC/` |
| `AOP_FASTQC_THREADS` | FastQC worker threads | `4` |
| `AOP_NULL_ITERATIONS` | Size-matched random draws per gene set and annotation family | `1000` |
| `AOP_NULL_SEED` | Reproducible matched-null seed | `20260825` |
| `AOP_NULL_DIRECTIONS` | Directions included in the matched-null analysis | `ALL` |
| `AOP_NETWORK_JSON` | AOPGraphExplorer JSON export used in the publication network figure | auto-discovered |
| `AOP_AOPXGENENET_RUN_ROOT` | Parent directory containing AOPxGeneNet/AOPxLINK runs | auto-discovered |
| `AOP_AOPXGENENET_RUN_DIR` | Explicit archived AOPxGeneNet/AOPxLINK run | auto-discovered |

For example, in a POSIX shell:

```bash
export AOP_CASE_STUDY_ID="study_b"
export AOP_COUNTS_FILE="data/study_b_counts.tsv"
export AOP_METADATA_FILE="data/study_b_metadata.csv"
export AOP_LITERATURE_FILE="data/study_b_literature.csv"
export AOP_CHEMICALS="CompoundA,CompoundB,PositiveControl"
export AOP_TIMEPOINT="48"
Rscript run_pipeline.R --steps=00,1,2,3,4,5,6,8,9,10,11,12
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
| 11 | Gene sets from 01/02, both branch-specific backgrounds, and AOPfingerprintR annotations |
| 12 | Available outputs from 01-11 plus optional AOPGraphExplorer/AOPxGeneNet exports |
| `fastqc` | Local FASTQ/FQ files, FastQC, and Jupyter; MultiQC is optional |

The current implementation is for human bulk RNA-seq: it uses human Ensembl
IDs, `org.Hs.eg.db`, human DoRothEA, KEGG `hsa`, Reactome `human`, and the human
AOP fingerprint resources. A non-human case study requires coordinated changes
to those annotation resources and organism arguments. The DESeq2 model is
`~ condition`; studies with batches, pairing, donor effects, or interactions
must extend the design formula in validation and differential-expression steps
rather than treating those variables as condition labels.

### Statistical families

- DESeq2 Wald-test P values are BH-adjusted across tested genes separately for
  each exposure-versus-control contrast.
- GO, KEGG, and Reactome over-representation results are BH-adjusted within each
  gene set and ontology/database collection. KE and AOP tests are likewise
  treated as separate annotation families.
- Step 11 uses one-sided hypergeometric tests, BH correction within each gene
  set × KE/AOP family, and an empirical P value based on the number of
  significant terms in 1,000 size-matched draws. Empirical P values are then
  BH-adjusted across evaluated sets within the KE and AOP families.
- WGCNA module-trait P values are nominal and exploratory. AOPxLINK
  module-AOP P values are BH-adjusted across the archived Ta4C3 test family.

Throughout the repository, `padj` denotes a BH-adjusted P value. “FDR < 0.05”
describes the corresponding decision threshold within the explicitly stated
test family; it is not a separate statistical test.

## Reproducibility and publication notes

- The random seed and shared paths are defined in `scripts/_config.R`.
- Raw data and generated outputs are not tracked. Document stable accessions,
  checksums, preprocessing, and download instructions for every released input.
- Review the chemical list and thresholds in `scripts/_config.R` before running
  a different dataset.
- The included GitHub Actions workflow checks R and Python syntax without
  downloading the full scientific software stack.
- [`CITATION.cff`](CITATION.cff), [`.zenodo.json`](.zenodo.json), and the
  [`Supporting Information manifest`](manuscript/SUPPORTING_INFORMATION_MANIFEST.md)
  provide release-ready citation and archive metadata. A Zenodo DOI exists only
  after an authorized GitHub release has actually been archived; no DOI is
  fabricated in this repository.

## License

No software license has been selected. Add a `LICENSE` file before public release
if you want others to be able to reuse or modify the code under stated terms.

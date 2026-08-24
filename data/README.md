# Input data

The pipeline expects input files under this directory, but data are ignored by
Git by default to prevent accidental redistribution of large or restricted
datasets. Review each source's license and consent requirements before changing
that policy.

## Required for differential expression

- `GSE200036_counts.tsv` (default): tab-separated raw integer counts, with human
  Ensembl gene IDs in the
  first column and sample names in the remaining columns.
- `GSE200036_metadata.csv` (default): sample metadata containing `sample_name` and
  `condition`. Sample names must match count-matrix columns. Control conditions
  may be named `CTRL`, `Control`, or `Medium`.

Different filenames can be supplied through `AOP_COUNTS_FILE` and
`AOP_METADATA_FILE`; they do not need to use the defaults above. Exposed
condition names must match the comma-separated `AOP_CHEMICALS` setting exactly.

## Required for the literature branch

- `Literature_DEGs.csv` (default): comma-separated table containing `Chemical` and a gene
  column such as `Gene`, `SYMBOL`, or `ENSEMBL`. `Direction` and `Organism`
  columns are recommended for provenance.
- `background_genes.xlsx`: one-column Ensembl background for the transcriptomic
  branch.
- `background_genes_GRCh38.xlsx`: one-column Ensembl background for the
  literature branch.

Optional step 8 requires `Direction` values that can be interpreted as UP/DOWN
or positive/negative numbers. It projects those literature signatures onto the
DESeq2 contrasts produced by step 1.

Override these paths with `AOP_LITERATURE_FILE`,
`AOP_BACKGROUND_GSE_FILE`, and `AOP_BACKGROUND_LIT_FILE`. Each background
workbook must contain a single human Ensembl-ID column and should represent the
genes that could have been detected in its branch, not only significant genes.

## Optional network inputs

Step 7 expects the following AOP-Wiki exports under `data/aop_wiki/`:

- `aop_ke_ec.tsv`: Key Event metadata.
- `aop_ke_ker.tsv`: Key Event relationships.
- `aop_ke_mie_ao.tsv`: Key Event roles and AOP membership.

Both step-07 implementations accept command-line path overrides; run either
script with `--help` for details.

# Supporting Information manifest

This manifest connects the manuscript claims to reproducible repository
artifacts. Generated `outputs/` are excluded from Git because they are large;
the archived release should deposit the listed outputs or an equivalent
checksum-verified results bundle.

## Curated evidence and provenance

- `MXene_data_template.xlsx`: screened-study metadata, gene evidence,
  literature DEG summaries, material characterization, and column definitions.
- `pubmed_chminsight_results_2025-09-06_19-25-56.xlsx`: unmodified PubMed
  ChemInsight export used for study screening.
- `notebooks/Fastq_QC.ipynb`: portable optional FastQC/MultiQC workflow.

## Supporting tables

- `outputs/Specificity_null/specificity_summary.csv`: observed and null KE/AOP
  counts, empirical P values, and empirical BH-adjusted P values.
- `outputs/Specificity_null/observed_term_level_enrichment.csv`: direct
  gene-to-KE/AOP sensitivity-analysis results.
- `outputs/Specificity_null/specificity_provenance.txt`: seed, iterations,
  backgrounds, correction families, and empirical statistic.
- `outputs/Publication_figures/publication_figure_manifest.csv`: figure source,
  dimensions, and generated-file inventory.

## Supporting figures

- `outputs/Publication_figures/FigureS1_AOP_specificity_null.png` and `.pdf`:
  observed enrichment versus 1,000 branch-specific, size-matched random gene
  sets.

## External application exports

- AOPGraphExplorer JSON export used for the focused/global AOP network panel.
- Archived AOPxGeneNet/AOPxLINK Ta4C3 run directory used for the exploratory
  module-trait and module-AOP panels.

## Release checklist

1. Run `Rscript run_pipeline.R --steps=00,1,2,3,4,5,6,8,9,10,11,12` in the
   restored `renv` environment.
2. Add the authorized AOPGraphExplorer and AOPxGeneNet/AOPxLINK exports to the
   archive bundle, preserving their original files.
3. Record SHA-256 checksums for inputs and generated result bundles.
4. Create a tagged GitHub release and enable Zenodo archiving.
5. Insert the resulting DOI into the manuscript and `CITATION.cff`; do not use a
   provisional or invented DOI.

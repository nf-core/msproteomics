# nf-core/msproteomics: Output

## Introduction

This document describes the output produced by the pipeline.
The directories listed below will be created in the results directory after the pipeline has finished.
All paths are relative to the top-level results directory specified with `--outdir`.

## Common Outputs

All workflow types produce the following:

### Pipeline Information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report with resource usage statistics.
  - `execution_timeline.html`: Timeline of task execution.
  - `execution_trace.txt`: Tab-separated trace of all tasks with resource metrics.
  - `pipeline_dag.dot` / `pipeline_dag.svg`: Pipeline DAG visualization.
  - `software_versions.yml`: Versions of all software used in the run.
  - `params.json`: Parameters used for the pipeline run.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides detailed execution reports for troubleshooting and reproducibility.

## DIA Workflow Output

DIA workflows (`--mode diann`) produce quantified protein and precursor matrices using DIA-NN.
DIA variants (e.g., phospho) are applied via `-c conf/variants/*.config` files.

### DIA-NN Results

<details markdown="1">
<summary>Output files</summary>

- `diann/`
  - `report.tsv`: Main DIA-NN output report with precursor-level quantification across all samples.
  - `report.stats.tsv`: Summary statistics (identifications, data completeness).
  - `report.pg_matrix.tsv`: Protein group-level quantification matrix.
  - `report.pr_matrix.tsv`: Precursor-level quantification matrix.
  - `report.gg_matrix.tsv`: Gene group-level quantification matrix.
  - `report.unique_genes_matrix.tsv`: Unique gene quantification matrix.
  - `*.speclib`: Spectral library generated or used by DIA-NN.

</details>

[DIA-NN](https://github.com/vdemichev/DiaNN) performs library-free or library-based DIA analysis with neural network-based scoring.

### MSstats Output

<details markdown="1">
<summary>Output files</summary>

- `msstats/`
  - `msstats_input.csv`: Formatted input for MSstats.
  - `msstats_results.csv`: Differential abundance analysis results (if contrasts are defined).

</details>

[MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html) provides statistical analysis for quantitative proteomics experiments.

### IQ (MaxLFQ) Quantification

<details markdown="1">
<summary>Output files</summary>

- `iq/`
  - `protein_maxlfq.tsv`: Protein-level MaxLFQ intensities.
  - `protein_maxlfq_log2.tsv`: Log2-transformed MaxLFQ intensities.

</details>

[IQ](https://github.com/tvpham/iq) implements the MaxLFQ algorithm for label-free protein quantification.

### MultiQC Report

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Aggregated QC report across all samples.
  - `multiqc_data/`: Parsed statistics from pipeline tools.

</details>

[MultiQC](http://multiqc.info) aggregates QC metrics from multiple tools into a single interactive HTML report.

### Interpreting DIA Results

Key metrics to check in `report.stats.tsv`:

- `Precursors.Identified`: typical range is 5,000-15,000 for a standard human DIA experiment, depending on instrument and method.
- `Proteins.Identified`: typically 3,000-8,000 for human plasma or cell lysate.
- High data completeness (>70% across samples) indicates good quantification quality.

The `report.pg_matrix.tsv` protein group matrix is the primary quantification output.
Columns represent samples; values are normalized quantities.
Missing values indicate the protein was not quantified in that sample.

For downstream analysis, import `report.pg_matrix.tsv` into tools like [Perseus](https://maxquant.net/perseus/), R, or Python for statistical testing, clustering, and visualization.

## DDA LFQ Workflow Output

The DDA LFQ workflow (`--mode fragpipe` without `--tmt_mode`) produces protein, peptide, and ion quantification using FragPipe tools.

### MSFragger Search Results

<details markdown="1">
<summary>Output files</summary>

- `msfragger/`
  - `*.pepXML`: Peptide-spectrum match results per sample.
  - `*.pin`: Percolator input files with rescoring features.

</details>

[MSFragger](https://msfragger.nesvilab.org/) performs ultrafast database search for peptide identification.

### Percolator Results

<details markdown="1">
<summary>Output files</summary>

- `percolator/`
  - `*.pout`: Percolator output with FDR-controlled PSMs per sample.

</details>

[Percolator](http://percolator.ms/) applies semi-supervised learning for PSM-level FDR control.

### Philosopher Filter and Report

<details markdown="1">
<summary>Output files</summary>

- Per-sample directories containing:
  - `psm.tsv`: Filtered PSM-level report.
  - `peptide.tsv`: Filtered peptide-level report.
  - `protein.tsv`: Filtered protein-level report.
  - `ion.tsv`: Filtered ion-level report.

</details>

[Philosopher](https://philosopher.nesvilab.org/) provides FDR filtering at protein, peptide, and PSM levels, then generates tabular reports.

### IonQuant Quantification

<details markdown="1">
<summary>Output files</summary>

- `ionquant/`
  - `combined_protein.tsv`: Combined protein-level quantification across all samples with MaxLFQ intensities.
  - `combined_peptide.tsv`: Combined peptide-level quantification.
  - `combined_ion.tsv`: Combined ion-level quantification.
  - `combined_modified_peptide.tsv`: Combined modified peptide quantification.
  - Per-sample `psm.tsv` files updated with quantification values.

</details>

[IonQuant](https://ionquant.nesvilab.org/) performs label-free quantification with match-between-runs and MaxLFQ normalization.

### Interpreting DDA LFQ Results

Key metrics to check:

- `combined_protein.tsv`: total protein count (typically 2,000-6,000 for a standard human experiment).
- Per-sample `protein.tsv` reports: compare protein counts across samples for consistency.
- `combined_ion.tsv`: ion-level data useful for more granular quantification analysis.

The `MaxLFQ Intensity` columns in `combined_protein.tsv` are the recommended quantification values.
These are normalized, label-free intensities suitable for cross-sample comparison.

Quality indicators:

- Consistent protein counts across samples (within 20% variance).
- Low percentage of single-peptide protein identifications.
- FDR columns should show values < 0.01 for high-confidence identifications.

For downstream analysis, use `combined_protein.tsv` with tools like [Perseus](https://maxquant.net/perseus/), [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html), or custom R/Python scripts.

## TMT Label Check Workflow Output

TMT label check workflows (`--mode fragpipe --tmt_mode labelcheck`) assess TMT labeling efficiency.

### TMT Labeling Efficiency Reports

<details markdown="1">
<summary>Output files</summary>

- `tmt_labelcheck/`
  - `tmt_labelcheck_report.html`: Interactive HTML report with labeling efficiency per sample and per channel.
  - `tmt_labelcheck_summary.tsv`: Tab-separated summary of labeling efficiency metrics.
  - `tmt_labelcheck_report.md`: Markdown summary report.

</details>

The TMT label check module analyzes PSM-level data to calculate TMT labeling efficiency, identifying under-labeled samples before committing to a full TMT quantification experiment.

### FragPipe Search Results

TMT label check workflows also produce the same MSFragger, Percolator, and Philosopher outputs described in the DDA LFQ section above, which are used as input to the labeling efficiency analysis.

### Interpreting TMT Label Check Results

The `tmt_labelcheck_report.html` is the primary QC output.
Key metrics:

- **Labeling efficiency**: percentage of PSMs with TMT modification detected. Values >95% indicate successful labeling; <90% suggests labeling issues.
- **Per-channel distribution**: all channels should show comparable PSM counts. Large imbalances may indicate mixing ratio problems.
- **Unlabeled peptide fraction**: high rates of unlabeled peptides may warrant re-labeling before proceeding to a full TMT quantification experiment.

If labeling efficiency is satisfactory, proceed to full TMT quantification with `--tmt_mode quant`.

## TMT Quantification Workflow Output

TMT quantification workflows (`--mode fragpipe --tmt_mode quant`) produce isobaric quantification results using TMTIntegrator.
Output includes the same MSFragger, Percolator, and Philosopher results described in the DDA LFQ section, plus TMTIntegrator quantification tables.

### Interpreting TMT Quantification Results

TMTIntegrator produces normalized, batch-corrected isobaric quantification.
The output tables contain ratio or intensity values per TMT channel per protein.
Use a reference channel (if included in your design) for normalization.

## Generic FragPipe Workflow Output

The generic FragPipe workflow produces output determined by the `.workflow` configuration file provided.
Output may include any combination of:

- MSFragger search results
- Percolator/PeptideProphet validation results
- ProteinProphet protein inference
- Philosopher filter and report files
- IonQuant or TMTIntegrator quantification
- PTMShepherd PTM analysis
- Spectral library files
- MSstats-compatible output

The exact output structure depends on which tools are enabled in the `.workflow` file.
Refer to the [FragPipe documentation](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) for details on each output type.

## Downstream Analysis

After pipeline completion, results can be analyzed with:

| Tool                                                                            | Use Case                                        | Input File                                                |
| ------------------------------------------------------------------------------- | ----------------------------------------------- | --------------------------------------------------------- |
| [Perseus](https://maxquant.net/perseus/)                                        | Statistical analysis, clustering, visualization | `combined_protein.tsv` or `report.pg_matrix.tsv`          |
| [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html) | Differential abundance testing                  | `msstats_input.csv` (DIA) or `combined_protein.tsv` (DDA) |
| [R/Bioconductor](https://www.bioconductor.org/)                                 | Custom statistical analysis                     | Any TSV output                                            |
| [Fragpipe-Analyst](https://fragpipe-analyst.nesvilab.org/)                      | Interactive visualization of FragPipe results   | FragPipe output directory                                 |

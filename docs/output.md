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
  - `execution_report_<timestamp>.html`: Nextflow execution report with resource usage statistics.
  - `execution_timeline_<timestamp>.html`: Timeline of task execution.
  - `execution_trace_<timestamp>.txt`: Tab-separated trace of all tasks with resource metrics.
  - `pipeline_dag_<timestamp>.html`: Pipeline DAG visualization.
  - `nf_core_msproteomics_software_mqc_versions.yml`: Versions of all software used in the run.
  - `params_<timestamp>.json`: Parameters used for the pipeline run.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides detailed execution reports for troubleshooting and reproducibility.

### MultiQC Report

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Software versions, parameter summary and methods description of the run.
  - `multiqc_data/`: Parsed statistics from pipeline tools.

</details>

[MultiQC](http://multiqc.info) aggregates QC metrics from multiple tools into a single interactive HTML report.

## DIA Workflow Output

DIA workflows (`--mode diann`) produce quantified protein and precursor matrices using DIA-NN.
DIA variants (e.g., phospho) are applied via `-c conf/variants/*.config` files.

### DIA-NN Results

<details markdown="1">
<summary>Output files</summary>

- `diann/`
  - `preliminary_analysis/`, `empirical_library/`, `individual_analysis/`: Outputs of the DIA-NN steps (spectral libraries, per-file `.quant` files, logs).
  - `diann_quantification/diann_report.tsv` (and `.parquet`): Main DIA-NN report with precursor-level quantification across all samples.
  - `diann_quantification/diann_report.stats.tsv`: Summary statistics (identifications, data completeness).
  - `diann_quantification/diann_report.pg_matrix.tsv`, `.pr_matrix.tsv`, `.gg_matrix.tsv`, `.unique_genes_matrix.tsv`: Protein group, precursor, gene group and unique gene matrices.
- `quantmsutils/`
  - `mztab/`: mzTab, MSstats and Triqler inputs converted from the DIA-NN report.
  - `mzml_statistics/`: Per-file mzML statistics (parquet).
- `pmultiqc/`: Proteomics QC report (`--enable_pmultiqc`).

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
  - `maxlfq.tsv`: Protein-level MaxLFQ intensities.
  - `contaminants_removed_*.tsv`: DIA-NN report without contaminant proteins (when `--contaminant_pattern` is set).

</details>

[IQ](https://github.com/tvpham/iq) implements the MaxLFQ algorithm for label-free protein quantification.

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

Pipeline mode (`--fragpipe_mode pipeline`) and the TMT label check also write FragPipe's run-folder files to `fragpipe/`, as FragPipe does in its work directory:

- `fragpipe-files.fp-manifest`: The LC-MS files with experiment, bioreplicate and data type.
- `fragpipe.workflow`: FragPipe's saved UI state for the run (every UI key with its default, run values and the tool-version header), written by FragPipe's own load and save code (`assets/fragpipe_ui_state/FragpipeUiState`, `FRAGPIPESAVEWORKFLOW`).
- `sdrf.tsv`: FragPipe's SDRF template (`workflow.misc.save-sdrf=true`, the default), written from the MSFragger keys of `fragpipe.workflow`. The fragment tolerance is FragPipe's first `New fragment_mass_tolerance` in the MSFragger output, read in manifest order as in FragPipe's per-file search.
- With `tab-run.write_sub_mzml=true` (FragPipe's first pass of a two-pass search):
  - `<run>_sub.mzML`: the scans of each LC-MS file that no PSM above `tab-run.sub_mzml_prob_threshold` explains, with MSFragger's calibrated peaks (`WRITESUBMZML`, FragPipe's WriteSubMzml).
  - `fragpipe-second-pass.workflow`: the run's `fragpipe.workflow` rewritten for the second pass (Percolator on, no mass calibration, MSFragger's optimized values from its output, FragPipe's `databases` folder, empty decoy tag and workdir), as FragPipe's finalizer writes it.
  - `fragpipe-files-second-pass.fp-manifest`: the manifest of the `fragpipe/<run>_sub.mzML` files.
  - `<group>/<run>_percolator.weights`: Percolator's feature weights, where the second pass finds them next to `<run>_sub.mzML` (FragPipe keeps them in its work dir's group folder).

Optional post-run steps that FragPipe runs from its Downstream tab, not as part of its run, are written under `downstream/`:

- `downstream/saintexpress/`: SAINTexpress (`saintexpress.run-saint-express=true`), one `saintexpress-<mode>/` per mode (`spc`, `int`) with the `inter/bait/prey.dat` inputs and SAINTexpress's `list.txt`, or `skipped_mode.txt` with FragPipe's error when the mode's `reprint.<mode>.tsv` is missing.

The exact output structure depends on which tools are enabled in the `.workflow` file.
Refer to the [FragPipe documentation](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) for details on each output type.

### Spectral Library (SpecLibGen)

<details markdown="1">
<summary>Output files</summary>

- `speclib/`
  - `easypqp_rt_reference_run.tsv`: EasyPQP reference-run retention times (automatic RT calibration).
  - `library.tsv`: Spectral library written by FragPipe's `gen_con_spec_lib.py` (EasyPQP library post-processed by `easypqp_lib_export('Spectronaut')`); this is the DIA-NN input.
  - `easypqp_files/`: `cmds.txt` with the FragPipe-SpecLib commands and the RT/IM alignment PDFs.
  - `filelist_speclibgen.txt`, `filelist_easypqp_library.txt`: Spectra and EasyPQP input lists.
  - `speclibgen_argv.txt`: The `gen_con_spec_lib.py` command, one argument per line.
  - `irt.tsv` / `im.tsv`, `*_temp-psm.tsv`: Calibration references and per-run PSM splits; with `speclibgen.keep-intermediate-files=true` also `*.psmpkl`, `*.peakpkl`, `*_run_peaks.tsv` and `easypqp_lib_openswath.tsv`.

</details>

### DIA-NN Quantification (FragPipe mode)

<details markdown="1">
<summary>Output files</summary>

- `dia_quant/`
  - `dia-quant-output/`: FragPipe's DIA-NN output directory: `report.tsv` (with propagated protein and gene information), `report.pr_matrix.tsv`, `report.pg_matrix.tsv`, `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv`, `report.stats.tsv`, `report-tsv.speclib` (DIA-NN's library, used by Skyline), `msstats.csv` (with `diann.generate-msstats=true`) and the FragReporter tables of the levels set in the workflow.
  - `experiment_annotation.tsv`: FragPipe's experiment annotation of the DIA runs.
  - `filelist_diann.txt`: DIA-NN inputs, largest first.
  - `diann_argv.txt`: The DIA-NN command, one argument per line.

</details>

## Downstream Analysis

After pipeline completion, results can be analyzed with:

| Tool                                                                            | Use Case                                        | Input File                                                |
| ------------------------------------------------------------------------------- | ----------------------------------------------- | --------------------------------------------------------- |
| [Perseus](https://maxquant.net/perseus/)                                        | Statistical analysis, clustering, visualization | `combined_protein.tsv` or `report.pg_matrix.tsv`          |
| [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html) | Differential abundance testing                  | `msstats_input.csv` (DIA) or `combined_protein.tsv` (DDA) |
| [R/Bioconductor](https://www.bioconductor.org/)                                 | Custom statistical analysis                     | Any TSV output                                            |
| [Fragpipe-Analyst](https://fragpipe-analyst.nesvilab.org/)                      | Interactive visualization of FragPipe results   | FragPipe output directory                                 |

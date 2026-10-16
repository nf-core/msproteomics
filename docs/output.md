# nf-core/msproteomics: Output

## Introduction

This document describes the output produced by the pipeline.
The directories listed below are created in the results directory after the pipeline has finished.
All paths are relative to the top-level results directory specified with `--outdir`.
Files are published with `--publish_dir_mode` (default `copy`).

In the FragPipe workflows the results directory takes the place of FragPipe's work directory: wherever FragPipe writes a file to its work directory root, the pipeline publishes it at the results root.
Placeholders used below:

- `<sample>`: the name of an LC-MS file without its extension.
- `<group>`: a FragPipe experiment/bioreplicate group (FragPipe's group folder in its work directory); the files of the blank group `""` are published at the top of the step directory, as FragPipe writes them to its work directory root.
- `fragpipe`: the name of the run-level tasks in FragPipe mode (for example `ionquant/fragpipe/`).

## Pipeline overview

Outputs common to all modes:

- [Samplesheet SDRF](#samplesheet-sdrf) - the samplesheet as an SDRF file
- [qpx export](#qpx-export) - DIA-NN, FragPipe LFQ and TMT-Integrator results as `.h5mu`, `.h5ad` and QPX Parquet
- [MultiQC](#multiqc) - software versions, parameters and methods description
- [Pipeline information](#pipeline-information) - report metrics generated during the workflow execution

DIA workflow (`--mode diann`):

- [File preparation](#file-preparation) - decompression logs, converted, indexed and peak-picked spectra
- [DIA-NN](#dia-nn) - in-silico library, preliminary analysis, empirical library and final quantification
- [quantms-utils](#quantms-utils) - DIA-NN configuration, mzML statistics and mzTab conversion
- [MSstats](#msstats) - MSstats LFQ analysis
- [IQ](#iq) - MaxLFQ protein quantification
- [pmultiqc](#pmultiqc) - proteomics QC report
- [Reference files](#reference-files) - DIA-NN in-silico library (`--save_reference`)

FragPipe workflows (`--mode fragpipe`, including the TMT label check):

- [FragPipe run files](#fragpipe-run-files) - manifest, saved workflow, SDRF template and write-sub-mzML files
- [FragPipe database](#fragpipe-database) - FragPipe's database update (`--save_reference`)
- [DIA-Umpire and diaTracer](#dia-umpire-and-diatracer) - DIA pseudo-spectra
- [MSFragger](#msfragger) - database search
- [Crystal-C](#crystal-c) - artifact removal
- [MSBooster](#msbooster) - deep-learning rescoring features
- [PSM validation](#psm-validation) - Percolator, PeptideProphet and PTMProphet
- [Protein inference](#protein-inference) - ProteinProphet, iProphet and Abacus
- [Philosopher filter and report](#philosopher-filter-and-report) - FDR filtering and the PSM, ion, peptide and protein tables
- [PTM-Shepherd and glyco analysis](#ptm-shepherd-and-glyco-analysis) - PTM-Shepherd, PairScans, O-Pair and MBG
- [IonQuant](#ionquant) - label-free, MS1 and isobaric quantification
- [TMT-Integrator](#tmt-integrator) - TMT reports
- [Spectral library and DIA-NN](#spectral-library-and-dia-nn) - EasyPQP library and DIA-NN quantification
- [FPOP](#fpop) - FPOP quantification and coADAPTr
- [Metaproteomics](#metaproteomics) - FragPipe's metaproteomics module
- [SAINTexpress](#saintexpress) - interaction scoring from FragPipe's Downstream tab
- [FragPipe headless](#fragpipe-headless) - all-in-one FragPipe run (`--fragpipe_mode headless`)
- [TMT label check](#tmt-label-check) - labeling efficiency reports

Each FragPipe step is published only when the `.workflow` file runs that tool.

## Common outputs

### Samplesheet SDRF

<details markdown="1">
<summary>Output files</summary>

- `sdrf/`
  - `<samplesheet>.sdrf.tsv`: the input samplesheet converted to SDRF, with the organism from `--organism`, the iProteome entry or the searched FASTA.

</details>

### qpx export

<details markdown="1">
<summary>Output files</summary>

- `qpx/`
  - `diann.h5mu`, `diann.proteins.h5ad`, `diann.*.parquet`, `diann.summary.tsv`: the DIA-NN final quantification (`--mode diann`).
  - `fragpipe.h5mu`, `fragpipe.proteins.h5ad`, `fragpipe.*.parquet`, `fragpipe.summary.tsv`: label-free IonQuant results (`--mode fragpipe` with IonQuant and without TMT-Integrator).
  - `fragpipe.tmt_<norm>.h5mu`, `fragpipe.tmt_<level>_<norm>.h5ad`, `fragpipe.tmt.summary.tsv`: TMT-Integrator abundance reports, one `.h5mu` per normalization with one modality per level (gene, protein, peptide, single-site, multi-site) and one `.h5ad` per level and normalization.

</details>

[qpx](https://github.com/bigbio/qpx) 1.1.4 (`qpxc convert diann|fragpipe`) writes the QPX Parquet dataset, and `bin/qpx_export.py` / `bin/tmt_to_anndata.py` build the [MuData](https://mudata.readthedocs.io/) and [AnnData](https://anndata.readthedocs.io/) files from it.
Observations are runs (DIA-NN, FragPipe LFQ; fractions of one experiment are merged after an identity check) or TMT samples, with the SDRF and annotation columns as `obs` metadata.
The `.h5mu` holds the `precursors` and `proteins` modalities; the `.h5ad` is the `proteins` modality.
DIA-NN: `proteins.X` = `PG.MaxLFQ` (pg_matrix), `precursors.X` = `Precursor.Normalised`, with `PG.Quantity` / `Precursor.Quantity` in `layers['raw']`.
FragPipe: `proteins.X` = `<experiment> MaxLFQ Intensity`, `precursors.X` = `<experiment> Intensity` of combined_ion.tsv.
Missing values are NaN, never 0, in a dense float32 matrix; `uns` records the quantity source and the missing-value convention, and `*.summary.tsv` the counts.
The export is on by default and a failed conversion fails the run; `--skip_qpx_export` turns it off.
A FragPipe LFQ export needs at least 3 experiment+bioreplicate combinations (IonQuant writes MaxLFQ intensities only then), so a run with fewer stops before the search; the TMT label check has no export.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: the data MultiQC parsed (software versions, parameter summary and methods description).
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project.
The report holds the software versions, the parameter summary and the methods description of the run.
No tool results feed it: in FragPipe mode (DDA LFQ, TMT and TMT Label Check) MultiQC gets no FragPipe inputs, only the parameter summary, the software versions and the methods description.

## DIA workflow

### File preparation

<details markdown="1">
<summary>Output files</summary>

- `decompress/`
  - `*_decompression.log`: log of the decompression of `.gz`, `.tar`, `.tar.gz` and `.zip` inputs (the decompressed spectra are copies of the inputs and are not published).
- `mzml/`
  - `<sample>.mzML`: spectra converted from Thermo `.raw` files by ThermoRawFileParser.
  - `out/<sample>.mzML`: re-indexed mzML files (`--reindex_mzml`).
  - `*_mzmlindexing.log`: log of the re-indexing.
- `tdf2mzml/`
  - `<sample>.mzML`, `*_conversion.log`: Bruker `.d` files converted to mzML (`--convert_dotd`).
- `openms/`
  - `picked/<sample>.mzML`: peak-picked (centroided) mzML files from OpenMS PeakPickerHiRes (`--openms_peakpicking`).
  - `*_pp.log`: log of the peak picking.

</details>

### DIA-NN

<details markdown="1">
<summary>Output files</summary>

- `diann/`
  - `preliminary_analysis/`: per-file DIA-NN analysis (`.quant` files, logs) for the mass accuracy and scan window calibration.
  - `empirical_library/`: the empirical spectral library assembled from the preliminary analysis (`empirical_library.parquet` or `.speclib`, logs).
  - `individual_analysis/`: per-file analysis against the empirical library (`.quant` files, logs).
  - `diann_quantification/`: the final quantification across all runs:
    - `diann_report.tsv` / `diann_report.parquet`: the main precursor-level report.
    - `diann_report.pr_matrix.tsv`, `diann_report.pg_matrix.tsv`, `diann_report.gg_matrix.tsv`, `diann_report.unique_genes_matrix.tsv`: precursor, protein group, gene group and unique gene matrices.
    - `diann_report.stats.tsv`: per-run identification statistics.
    - `diann_report.log.txt` and the other DIA-NN side files (`.manifest.txt`, `.protein_description.tsv`, `.speclib`).

</details>

[DIA-NN](https://github.com/vdemichev/DiaNN) performs library-free or library-based DIA analysis with neural network-based scoring.
DIA variants (for example phospho) are applied with `-c conf/variants/*.config`.

### quantms-utils

<details markdown="1">
<summary>Output files</summary>

- `quantmsutils/`
  - `config/`: the DIA-NN configuration generated from the SDRF (`diann_config.cfg`) and its log.
  - `mzml_statistics/`: per-file mzML statistics (`<sample>_ms_info.parquet`, `<sample>_ms2_info.parquet`, `<sample>_feature_info.parquet`, logs).
  - `mztab/`: mzTab, MSstats and Triqler inputs converted from the DIA-NN report (`*.mzTab`, `*_msstats_in.csv`, `*_triqler_in.tsv`, `experimental_design_*`).

</details>

### MSstats

<details markdown="1">
<summary>Output files</summary>

- `msstats/`
  - `*.csv`: MSstats results (protein-level summaries and the comparisons when contrasts are defined).
  - `*.pdf`: MSstats plots (when they could be generated).
  - `*.log`: MSstats log.

</details>

[MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html) provides statistical analysis for quantitative proteomics experiments.
It is skipped with `--skip_post_msstats`.

### IQ

<details markdown="1">
<summary>Output files</summary>

- `iq/`
  - `maxlfq.tsv`: protein-level MaxLFQ intensities.
  - `contaminants_removed_*.tsv`: the DIA-NN report without contaminant proteins (when `--contaminant_pattern` is set).

</details>

[IQ](https://github.com/tvpham/iq) implements the MaxLFQ algorithm for label-free protein quantification.

### pmultiqc

<details markdown="1">
<summary>Output files</summary>

- `pmultiqc/`
  - `*.html`: proteomics QC report.
  - `*_data/`: data of the report.
  - `*.db`: report database (optional).

</details>

[pmultiqc](https://github.com/bigbio/pmultiqc) is published with `--enable_pmultiqc`.

### Reference files

<details markdown="1">
<summary>Output files</summary>

- `proteome/diann/`: the DIA-NN in-silico spectral library predicted from the FASTA, with `--save_reference`.

</details>

## FragPipe workflows

### FragPipe run files

<details markdown="1">
<summary>Output files</summary>

- `fragpipe/`
  - `fragpipe-files.fp-manifest`: the LC-MS files with experiment, bioreplicate and data type (`--fragpipe_mode pipeline`; not written by the TMT label check).
  - `fragpipe.workflow`: FragPipe's saved workflow of the run (every setting with its value and the tool-version header), written by FragPipe's own code.
  - `sdrf.tsv`: FragPipe's SDRF template, written when `workflow.misc.save-sdrf=true` (the default).
  - With `tab-run.write_sub_mzml=true` (the first pass of FragPipe's two-pass search):
    - `<sample>_sub.mzML`: the scans of each LC-MS file that no PSM above `tab-run.sub_mzml_prob_threshold` explains.
    - `fragpipe-second-pass.workflow`: the run's workflow rewritten for the second pass.
    - `fragpipe-files-second-pass.fp-manifest`: the manifest of the `_sub.mzML` files.
    - `<group>/<sample>_percolator.weights`: Percolator's feature weights, where the second pass finds them.
- `experiment_annotation.tsv`: FragPipe's experiment annotation of the run (IonQuant or, for isobaric labels, LabelQuant), at the results root.

</details>

### FragPipe database

<details markdown="1">
<summary>Output files</summary>

- `proteome/fragpipe/`
  - `<database>_philosopher.fasta`: the database after FragPipe's database update (decoys with `--add_philosopher_decoy`, contaminants with `--add_philosopher_contaminants`).
  - `philosopher_database.log`: log of `philosopher database`.

</details>

Published with `--save_reference`.

### DIA-Umpire and diaTracer

<details markdown="1">
<summary>Output files</summary>

- `diaumpire/<sample>/`
  - `<sample>_Q1.mzML`, `<sample>_Q2.mzML`, `<sample>_Q3.mzML` (or `.mgf`): DIA-Umpire pseudo-MS/MS spectra by quality tier.
  - `umpire-se.params`, `diaumpire.log`: DIA-Umpire parameters and log.
  - `diaumpire_se.log`: DIA-Umpire SE's own log, which FragPipe keeps in its work dir.
- `diatracer/<sample>/`
  - `<sample>_diatracer.mzML`: diaTracer pseudo-MS/MS spectra of diaPASEF data.
  - `diatracer.log`: diaTracer log.

</details>

### MSFragger

<details markdown="1">
<summary>Output files</summary>

- `msfragger/`
  - `fragger.params` (and `fragger_dda.params`, `fragger_dia.params`, `fragger_gpfdia.params`, `fragger_dda_plus.params` for the data types the run searches): FragPipe's MSFragger parameter files, with the run's database.
  - `<sample>/`: per-file search results:
    - `<sample>.pepXML`: peptide-spectrum matches.
    - `<sample>.pin`: Percolator input with the search features.
    - `<sample>_scores_histogram.tsv`, `fragger.params`, `msfragger.log`: score histogram, parameters and log of the search.
  - `<sample>_uncalibrated.mzML`, `<sample>_calibrated.mzML`, `<sample>_uncalibrated.mzBIN`, `<sample>.mzBIN_calibrated`: the spectra files MSFragger writes next to the LC-MS file (`.raw` or `.d` input, mass calibration, `write_calibrated_mzml`).
  - `<sample>.mzBIN_all`: all spectra for Write sub mzML (`write_mzbin_all = 1`).

</details>

[MSFragger](https://msfragger.nesvilab.org/) performs ultrafast database search for peptide identification.
In a split-database search (`msfragger.misc.slice-db` > 1) with mass calibration, `msfragger/` holds the `<sample>_uncalibrated.mzBIN` and, for `.raw` or `.d` input, the `<sample>_uncalibrated.mzML` the calibration search writes next to the LC-MS file; the calibrated spectra and parameters are not published, as FragPipe's `msfragger_pep_split.py` keeps them in a temporary directory that it deletes at the end.
With a split database, the merged results of the database chunks are published in `msfragger/<sample>/` like a single search.

### Crystal-C

<details markdown="1">
<summary>Output files</summary>

- `crystalc/<sample>/`
  - `*_c.pepXML`: the search results after artifact removal.
  - `crystalc.log`: Crystal-C log.

</details>

### MSBooster

<details markdown="1">
<summary>Output files</summary>

- `msbooster/`
  - `msbooster/`:
    - `<sample>_edited.pin`: Percolator input with MSBooster's deep-learning features.
    - `msbooster_params.txt`, `msbooster_msbooster.log`: MSBooster parameters and log.
  - `MSBooster/` (or `<group>/MSBooster/`): MSBooster's predictions and plots (`spectraRT.tsv`, `spectraRT_full.tsv`, `MSBooster_plots/`).

</details>

### PSM validation

<details markdown="1">
<summary>Output files</summary>

- `percolator/<sample>/`
  - `interact-<sample>.pep.xml`: Percolator-validated PSMs as pepXML.
  - `<sample>_target_psms.tsv`, `<sample>_decoy_psms.tsv`: Percolator target and decoy PSMs (optional).
  - `percolator.log`: Percolator log.
- `peptideprophet/<group>/`
  - `interact*.pep.xml`: PeptideProphet-validated PSMs.
  - `peptideprophet.log`: PeptideProphet log.
- `ptmprophet/<sample>/`
  - `*.mod.pep.xml`: PSMs with PTMProphet site localization.
  - `ptmprophet.log`: PTMProphet log.

</details>

[Percolator](http://percolator.ms/) applies semi-supervised learning for PSM-level FDR control.
PeptideProphet and PTMProphet run through [Philosopher](https://philosopher.nesvilab.org/).

### Protein inference

<details markdown="1">
<summary>Output files</summary>

- `proteinprophet/fragpipe/`
  - `combined.prot.xml`: ProteinProphet protein inference over all groups.
  - `filelist_proteinprophet.txt`, `proteinprophet.log`: input list and log.
- `iprophet/fragpipe/`
  - `combined.pep.xml`, `iprophet.log`: iProphet results of a multi-experiment report.
- `abacus/fragpipe/`
  - `combined_protein.tsv`, `combined_peptide.tsv`: Abacus multi-experiment protein and peptide tables.
  - `reprint.spc.tsv`, `reprint.int.tsv`: REPRINT spectral-count and intensity tables.

</details>

### Philosopher filter and report

<details markdown="1">
<summary>Output files</summary>

- `philosopher/`
  - `<group>/` (the blank group's files at the top of `philosopher/`):
    - `psm.tsv`: PSM table.
    - `ion.tsv`: ion table.
    - `peptide.tsv`: peptide table.
    - `protein.tsv`: protein table.
    - `philosopher_filter.log`, `filter.log`: logs of `philosopher filter` and `philosopher report`.
  - `export_matched_fragments.log`: log of FragPipe's "Export matched fragments" (when enabled).

</details>

[Philosopher](https://philosopher.nesvilab.org/) filters the PSMs, peptides and proteins at the workflow's FDR and writes the tables.
When FreeQuant or LabelQuant run, the tables come from the report that follows them.
FreeQuant and LabelQuant only update Philosopher's `.meta` workspace, so they publish nothing of their own.
As in FragPipe, `philosopher workspace --clean` removes each group's `.meta` workspace after its last Philosopher command (`PHILOSOPHERWORKSPACECLEAN`, or `PHILOSOPHERREPORT` after FreeQuant or LabelQuant): no published directory holds one.
With "Export matched fragments" (or FPOP coADAPTr), `psm.tsv` is the annotated (or edited) version, as FragPipe rewrites it in place.

### PTM-Shepherd and glyco analysis

<details markdown="1">
<summary>Output files</summary>

- `ptmshepherd/`
  - `ptm-shepherd-output/`: PTM-Shepherd's results (`global.modsummary.tsv`, `global.profile.tsv`, `glyco_masses_list.txt`, ...).
  - `shepherd.config`, `ptmshepherd.log`: PTM-Shepherd configuration and log.
  - `groups/<group>/`: the group's tables with the `psm.tsv` PTM-Shepherd rewrote.
- `pairscans/`
  - `<sample>.pairs`, `<sample>.pairscans.log`: paired scans for O-Pair.
- `opair/<group>/`
  - the group's tables with the `psm.tsv` O-Pair rewrote, `filelist_opair.txt` and `opair.log`.
- `mbg/`
  - `<group>/`: the group's tables with the `psm.tsv` MBG rewrote, and `mbg.log`.
  - `fragpipe-files.fp-manifest`: the manifest MBG reads.

</details>

### IonQuant

<details markdown="1">
<summary>Output files</summary>

- `ionquant/fragpipe/`: label-free quantification:
  - `combined_protein.tsv`, `combined_peptide.tsv`, `combined_modified_peptide.tsv`, `combined_ion.tsv`: combined tables across all runs with MaxLFQ intensities.
  - `combined_site_*.tsv`: site reports (optional).
  - `<group>/`: the group's tables with the quantified `psm.tsv`.
  - `reprint.spc.tsv`, `reprint.int.tsv`: REPRINT tables.
  - `filelist_ionquant.txt`, `modmasses_ionquant.txt`, `ionquant.log`: IonQuant inputs and log.
- `ionquant_ms1/fragpipe/`: IonQuant MS1 quantification of a TMT run (same layout).
- `ionquant_isobaric/fragpipe/`: IonQuant isobaric quantification of a TMT run (same layout).

</details>

[IonQuant](https://ionquant.nesvilab.org/) performs label-free quantification with match-between-runs and MaxLFQ normalization, and isobaric reporter-ion quantification.
The `MaxLFQ Intensity` columns of `combined_protein.tsv` are the recommended label-free protein quantities.

### TMT-Integrator

<details markdown="1">
<summary>Output files</summary>

- `tmt-report/`
  - `abundance_*.tsv`, `ratio_*.tsv`: TMT-Integrator abundance and ratio reports by level (gene, protein, peptide, site) and normalization.
- `tmt-integrator-conf.yml`: TMT-Integrator configuration.
- `MSstatsTMT_annotation.csv`: MSstatsTMT annotation (optional).
- `tmt-report-unmod/`, `tmt-integrator-conf-unmod.yml`: the second TMT-Integrator run on unmodified peptides (FPOP with TMT).

</details>

These files are at the results root, as in FragPipe's work directory.

### Spectral library and DIA-NN

<details markdown="1">
<summary>Output files</summary>

- `speclib/`
  - `library.tsv`: the spectral library written by FragPipe's `gen_con_spec_lib.py` (EasyPQP); the DIA-NN input.
  - `easypqp_rt_reference_run.tsv`: EasyPQP reference-run retention times (automatic RT calibration).
  - `easypqp_im_reference_run.tsv`, `easypqp_im_alignment_<run>.alignment_pkl`: EasyPQP reference-run ion mobilities and the ion mobility alignment of each run (ion mobility data, e.g. diaPASEF).
  - `easypqp_files/`: `cmds.txt` with the EasyPQP commands and the RT/IM alignment PDFs.
  - `filelist_speclibgen.txt`, `filelist_easypqp_library.txt`, `speclibgen_argv.txt`: spectra and EasyPQP input lists and the `gen_con_spec_lib.py` command.
  - `irt.tsv` / `im.tsv`, `*_temp-psm.tsv`: calibration references and per-run PSM splits; with `speclibgen.keep-intermediate-files=true` also `*.psmpkl`, `*.peakpkl`, `*_run_peaks.tsv` and `easypqp_lib_openswath.tsv`.
- `dia_quant/`
  - `dia-quant-output/`: FragPipe's DIA-NN output: `report.tsv`, `report.pr_matrix.tsv`, `report.pg_matrix.tsv`, `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv`, `report.stats.tsv`, `report-tsv.speclib`, `msstats.csv` (with `diann.generate-msstats=true`) and the FragReporter tables of the levels set in the workflow.
  - `experiment_annotation.tsv`: FragPipe's experiment annotation of the DIA runs.
  - `filelist_diann.txt`, `diann_argv.txt`: the DIA-NN inputs (largest first) and the DIA-NN command.

</details>

### FPOP

<details markdown="1">
<summary>Output files</summary>

- `fpop/fragpipe/`
  - `FPOP_peptides.tsv`, `FPOP_sites.tsv`: FPOP peptide and site quantification.
- `fpop_coadaptr/`
  - `fpop_coadaptr.log`: FPOP coADAPTr log.

</details>

FPOP coADAPTr edits the final `psm.tsv` files (or DIA-NN's `report.tsv`) in place, as FragPipe does: the edited tables are published at the path of the step that wrote them (`ionquant/`, `ionquant_isobaric/`, `philosopher/` or `dia_quant/`).

### Metaproteomics

<details markdown="1">
<summary>Output files</summary>

- `metaproteomics/fragpipe/`
  - `refined_db.fasta`, `combined_PSM_species.tsv`: FragPipe's metaproteomics results.

</details>

### SAINTexpress

<details markdown="1">
<summary>Output files</summary>

- `downstream/saintexpress/saintexpress-<mode>/` (`spc`, `int`)
  - `inter.dat`, `bait.dat`, `prey.dat`: SAINTexpress inputs from the REPRINT tables.
  - `list.txt`, `saintexpress.log`: SAINTexpress results and log.
  - `skipped_mode.txt`: FragPipe's error message when the mode's `reprint.<mode>.tsv` is missing.

</details>

SAINTexpress runs with `saintexpress.run-saint-express=true`.
It is a step of FragPipe's Downstream tab, not of FragPipe's run, so it is published under `downstream/`.

### FragPipe headless

<details markdown="1">
<summary>Output files</summary>

- `fragpipe/`
  - `results/`: FragPipe's work directory of the headless run (`combined_protein.tsv`, `combined_peptide.tsv`, `combined_ion.tsv`, per-group folders, logs and every other file FragPipe writes).
  - `manifest.fp-manifest`: the manifest FragPipe headless ran with.
  - `raw_files/*annotation.txt`: the isobaric channel annotations next to the LC-MS files (TMT).

</details>

`--fragpipe_mode headless` runs FragPipe in one task instead of the per-tool steps above.

### TMT label check

<details markdown="1">
<summary>Output files</summary>

- `qc/`: labeling efficiency from the PSMs of every LC-MS file:
  - `report.html`: HTML report with the labeling efficiency per sample.
  - `report.md`: Markdown version of the report.
  - `labeling_summary.tsv`: labeling efficiency metrics with the pass/warn/fail status.
  - `per_sample_efficiency.csv`: labeling efficiency per sample.
- `qc_combined/`: the same files from IonQuant's `combined_modified_peptide.tsv`.

</details>

The TMT label check (`--mode fragpipe --tmt_mode labelcheck`) also publishes the FragPipe steps its workflow runs (`msfragger/`, `percolator/`, `proteinprophet/`, `philosopher/`, `ionquant/fragpipe/`, `fragpipe/`).
A labeling efficiency of 95% or more passes, 85% or more warns, and less fails.
If the labeling efficiency is satisfactory, proceed to TMT quantification with `--tmt_mode quant`.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and `pipeline_dag_<timestamp>.html`.
  - Parameters used by the pipeline run: `params_<timestamp>.json`.
  - Software versions of the run: `nf_core_msproteomics_software_mqc_versions.yml`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline.
This will allow you to troubleshoot errors with the running pipeline, and also provide you with other information such as launch commands, run times and resource usage.

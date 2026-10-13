# nf-core/msproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- FragPipe 24.0 parity for `--mode fragpipe` (DDA LFQ, TMT label check, TMT quantification, generic `.workflow` files): FragPipe's tool graph (FragpipeRun.java) as modular subworkflows (`fragpipe_database`, `fragpipe_convert`, `fragpipe_search`, `fragpipe_validate`, `fragpipe_inference`, `fragpipe_philosopher_quant`, `fragpipe_glyco`, `fragpipe_ptm`, `fragpipe_quant`, `fragpipe_speclib`, `fragpipe_export`), one local module per FragPipe tool, and every argv and config file built from the `.workflow` file by `PARSEFRAGPIPEWORKFLOW` and Groovy helpers (`fragpipe_utils*.nf`) exactly as FragPipe builds them. Includes split-database search with calibration, write-sub-mzML second pass, glyco (PTM-Shepherd, MBG, O-Pair), FPOP, metaproteomics, SAINTexpress and DIA-NN spectral-library quantification; `--fragpipe_mode headless` runs FragPipe itself.
- iProteomes (`--proteome`, `conf/iproteomes.config`): versioned reference proteomes with hosted DIA-NN libraries and FragPipe decoy databases.
- Parity tests: parser and stub tests against FragPipe 24.0 behaviour, and golden tests (tag `golden`) that compare module commands and outputs with captured FragPipe 24.0 runs on public datasets (PXD001468, PXD017646, PXD030670, PXD047839, CPTAC).

### Fixed

- FragPipe tools run from both FragPipe distributions: the MSFragger, IonQuant and diaTracer jars are found in the task as FragPipe does (TabConfig.getJarPath: `<tool>-<version>.jar` under the tools folder, with or without an edition word), and batmass-io and MSBooster are the jars the image's FragPipe build names (ToolingUtils.BATMASS_IO_JAR, CmdMSBooster.JAR_MSBOOSTER_NAME). A missing jar stops the task with an error naming the tools folder.
- `IQ`: no longer fails with "Missing output file(s) `maxlfq.tsv`" when no report row passes the FDR filter (`Q.Value`, `PG.Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`; e.g. `PG.Q.Value = 1` for every row on a 500-protein FASTA): `process_long_format()` wrote no file and R exited 0. A header-only `maxlfq.tsv` (`Protein.Group`, `Protein.Names`, `Genes`, the samples) is now written with a `WARNING: IQ: maxlfq.tsv: no report row passes ...` line giving the filter and the minimum q-values. A report with no rows, or with only contaminant rows, fails with `ERROR: IQ: ...`, and so does a run where rows pass but no `maxlfq.tsv` is written. New integration tests reuse the `DIANNR` 10-row report fixture; the `DIANNR` test now reads the task's `.command.err` (not the Nextflow console) for the warning.
- `DIANNR`: no longer fails with `subscript out of bounds` in `diann_maxlfq` when no report row passes the FDR filter (e.g. DIA-NN reports `PG.Q.Value = 1` for every row on a small search space such as a 500-protein FASTA). `diann_maxlfq` is only called when rows pass; otherwise, and for the `diann_matrix` outputs that were silently written as 0-byte files, a header-only table (the sample columns) is written with a `WARNING: DIANNR: <file>: no report row passes ...` line. A report with no rows, or with only contaminant rows, fails with `ERROR: DIANNR: ...`; an empty table from rows that did pass fails too. New integration tests with a 10-row report fixture.
- `MERGE_SPLIT_SEARCH`: the task fails at once with `ERROR: MERGE_SPLIT_SEARCH: no MSFragger*.jar found` (naming the searched `msfragger_dir` and `ext.fragpipe_tools_dir`) when no MSFragger jar is found; before, `java -jar` got an empty jar path and reported a misleading "Invalid or corrupt jarfile <histogram>.tsv". `merge_split_search.py` also rejects a `--msfragger_cmd` whose `-jar` path is missing.
- `MERGE_SPLIT_SEARCH` tests: the integration test runs in the FragPipe image given by `FRAGPIPE_CONTAINER` (it ran on the host before); a new stub-tagged test covers the missing-jar error; the stub snapshot no longer records the host Python version.
- `FRAGPIPE_HEADLESS`: `manifest_content` is used again and is authoritative (file, experiment, bioreplicate, data type), so allinone runs keep one experiment per sample, empty bioreplicates for fractions, and `DDA+`; before, every LFQ file fell back to `experiment1` / `1` / `DDA`.
- `FRAGPIPE_HEADLESS`: file names are matched exactly against the manifest (or `file_experiment_map`) first column; the substring lookup (`grep -F`) gave `A1.mzML` the experiment of `XA1.mzML`.
- `FRAGPIPE_HEADLESS`: the task fails loudly when a staged file has no row, a row names an unstaged file, a file is listed twice, a TMT file's experiment is not in `annotation_content`, or both `manifest_content` and `file_experiment_map` are empty (TMT files were silently dropped before).
- `FRAGPIPE_HEADLESS`: `versions.yml` no longer carries a stray `END_VERSIONS` line when a multi-line value was pasted into the script (inputs are now passed base64-encoded, so no value can break indentation or expand in bash).
- `nf-core pipelines lint`: `params.diann_additional_analysis_args` defaults to `null` (schema has no default); nf-core reads a config `''` as `null`, so `''` could never match the schema. All consumers use `?: ''`, so behaviour is unchanged.
- `nf-core pipelines lint`: removed the `max_cpus` / `max_memory` / `max_time` schema params and the `test` profile's `params.max_*`; resource caps are set by `process.resourceLimits`, and nothing read the params.
- `nf-core pipelines lint`: removed `params.proteomes = null` from `conf/reference_proteomes_ignored.config`; no code reads `params.proteomes`.
- `utils_nfcore_pipeline`: tests match the installed nf-core/modules commit (`tests/main.nf.test` restored; stale `tests/main.workflow.nf.test` and its snapshot removed).
- Reference proteomes: the documented UniProt download now works. `conf/reference_proteomes.config` defines `params.databases` (keyed by organism, e.g. `'Homo sapiens'`, with a `database` URL), which is what the workflows read; before it defined an unused `params.proteomes` and `params.proteomes_ignore` defaulted to `true`, so without `--database` every run failed. Mirrors nf-core iGenomes: `proteomes_ignore` defaults to `false`, `proteomes_base` holds the UniProt URL prefix, `conf/reference_proteomes_ignored.config` sets `params.databases = [:]`, and `databases` is in `validation.defaultIgnoreParams`. `--database` stays authoritative; an organism with no proteome and no `--database` fails with `No protein database for organism '<organism>': ...` listing the available organisms (`getProteinDatabase`, new function tests).
- Database channels use `Channel.fromPath(..., glob: false)`; with globbing, the `?` in a UniProt REST URL dropped the query string.

- Sync with nf-core/msproteomics `dev` (the FragPipe 24.0 parity rewrite renamed or replaced several of the files above; each fix is kept in the new code):
  - `IQ`: the no-row-passes-FDR handling is in `modules/local/iq/templates/iq.r`; its regression report moved to `modules/local/iq/tests/diann_report_pgq1.tsv` (`DIANNR` was removed).
  - `MERGESPLITSEARCH` (was `MERGE_SPLIT_SEARCH`): the tool paths are built in Groovy, so the missing-jar check is in `merge_split_search.py`; it runs before any merge work and takes the jar as the first non-option token after `-jar`, since FragPipe puts JVM options there (`msfragger_pep_split.py:73`).
  - `FRAGPIPEHEADLESS` (was `FRAGPIPE_HEADLESS`): the manifest is built from structured rows (`headlessLayout`), never by substring lookup; `checkHeadlessRows` stops the run on a staged file without a row, a row naming no staged file, a file listed twice, no rows, or a TMT group without annotation channels. `versions.yml` is written with `printf`.
  - `glob: false` on every database path (`resolveFragPipeDatabaseChannels`, the DIA and FragPipe workflows); reference proteomes come from iProteomes (`--proteome`), which replace `params.databases`.
  - `test_tmt` sets `tmt_type = 'TMT6'` (TMT label check needs it to pick its workflow).
  - `utils_nfcore_pipeline` is the unpatched nf-core/modules version (no local `imNotification`), so its patch file is removed.

### Changed

- `FRAGPIPE_HEADLESS`: single-plex TMT files are placed in `raw_files/<experiment>/` like multi-plex ones.
- `FRAGPIPE_HEADLESS` stub writes `results/fragpipe-files.fp-manifest` and `results/experiment_annotation.tsv`.
- FragPipe helper functions (`shouldRunTool`, `getToolArgs`, `getToolModmasses`, `getToolField`, `getToolReportArgs`, `generateFragpipeManifest`) moved from `subworkflows/local/fragpipe_utils.nf` into `subworkflows/local/utils_nfcore_msproteomics_pipeline/main.nf`; `nf-core pipelines lint` crashed (`IndexError` in the subworkflow `main_nf` check) on the functions-only file. The unused duplicate `lib/FragpipeUtils.nf` is removed.

## v1.0.0 - Initial Release

### Added

- DIA workflow (DIA-NN) with MSstats and MaxLFQ quantification
- DDA LFQ workflow (FragPipe: MSFragger, MSBooster, Percolator, IonQuant)
- TMT Label Check workflow (FragPipe: TMT labeling efficiency QC)
- Generic FragPipe workflow (configurable via .workflow files)
- Multi-instrument support (ASC, FLX)
- SDRF input format support
- Pre-built database selection for human, mouse, and yeast
- nf-core native pipeline structure

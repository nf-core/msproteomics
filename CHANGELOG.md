# nf-core/msproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed

- `FRAGPIPE_HEADLESS`: `manifest_content` is used again and is authoritative (file, experiment, bioreplicate, data type), so allinone runs keep one experiment per sample, empty bioreplicates for fractions, and `DDA+`; before, every LFQ file fell back to `experiment1` / `1` / `DDA`.
- `FRAGPIPE_HEADLESS`: file names are matched exactly against the manifest (or `file_experiment_map`) first column; the substring lookup (`grep -F`) gave `A1.mzML` the experiment of `XA1.mzML`.
- `FRAGPIPE_HEADLESS`: the task fails loudly when a staged file has no row, a row names an unstaged file, a file is listed twice, a TMT file's experiment is not in `annotation_content`, or both `manifest_content` and `file_experiment_map` are empty (TMT files were silently dropped before).
- `FRAGPIPE_HEADLESS`: `versions.yml` no longer carries a stray `END_VERSIONS` line when a multi-line value was pasted into the script (inputs are now passed base64-encoded, so no value can break indentation or expand in bash).
- `nf-core pipelines lint`: `params.diann_additional_analysis_args` defaults to `null` (schema has no default); nf-core reads a config `''` as `null`, so `''` could never match the schema. All consumers use `?: ''`, so behaviour is unchanged.
- `nf-core pipelines lint`: removed the `max_cpus` / `max_memory` / `max_time` schema params and the `test` profile's `params.max_*`; resource caps are set by `process.resourceLimits`, and nothing read the params.
- `nf-core pipelines lint`: removed `params.proteomes = null` from `conf/reference_proteomes_ignored.config`; no code reads `params.proteomes`.
- `utils_nfcore_pipeline`: tests match the installed nf-core/modules commit (`tests/main.nf.test` restored; stale `tests/main.workflow.nf.test` and its snapshot removed).

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

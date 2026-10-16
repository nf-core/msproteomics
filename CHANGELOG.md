# nf-core/msproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [2026-10-14]

### Added

- `main.nf` follows the nf-core template: `PIPELINE_INITIALISATION` (nf-schema parameter and samplesheet validation, `--help`, `--version`), the `NFCORE_MSPROTEOMICS` named workflow, a pipeline-level `MULTIQC` with software versions from the `versions` topic, and `PIPELINE_COMPLETION` (summary and email); none of these ran before.
- `PMULTIQC` runs in the DIA workflow when `--enable_pmultiqc` is set (the default); the module existed but was never called.
- iProteomes: `--proteome <key>` selects a versioned reference proteome from `conf/iproteomes.config` (paths under `--iproteomes_base`) with its FASTA and, when the run's settings match, a hosted DIA-NN library and FragPipe decoy database.
- `--save_reference` publishes the DIA-NN library and FragPipe database built by the pipeline to `<outdir>/proteome/`.
- `--add_philosopher_decoy` and `--add_philosopher_contaminants` run FragPipe's database update (`philosopher database`) on `--database` before a FragPipe search.
- Samplesheet columns `experiment`, `bioreplicate` and `data_type` set FragPipe's experiment, bioreplicate and data type per LC-MS file.
- The samplesheet `spectra` column accepts `.mzXML`, `.mgf` and `.mzBIN` files and `.gz`, `.tar`, `.tar.gz` and `.zip` archives (`DECOMPRESS` unpacks them in the DIA workflow).
- `--tmt_type TMT35`.
- FragPipe tools follow FragPipe 24.0's tool graph: `CHECKCENTROID`, `CHECKDBCONFIG`, iProphet (`PHILOSOPHERIPROPHET`), Abacus (`PHILOSOPHERABACUS`), matched-fragment export (`EXPORTMATCHEDFRAGMENTS`), write-sub-mzML second pass (`WRITESUBMZML`, `FRAGPIPESECONDPASS`), split-database calibration (`SORTFASTA`, `MSFRAGGERAPPLYCALIBRATION`), O-Pair scan pairing (`PAIRSCANS`), SAINTexpress input from IonQuant (`REPRINTTOSAINTINPUT`), FPOP coADAPTr (`FPOPCOADAPTR`) and DIA-NN quantification of a FragPipe spectral library (`DIANNFRAGPIPE`).
- New local subworkflows `fragpipe_philosopher_quant` and `fragpipe_matched_fragments`.
- `--mode fragpipe` publishes `fragpipe/fragpipe-files.fp-manifest`, `fragpipe/fragpipe.workflow` (`FRAGPIPESAVEWORKFLOW`) and, with `workflow.misc.save-sdrf=true`, `fragpipe/sdrf.tsv` (`FRAGPIPESDRF`).
- `--fragpipe_threads` sets the FragPipe thread count; `--metaproteomics_names_dmp` and `--metaproteomics_nodes_dmp` give the NCBI taxonomy files for FP-Meta.
- `test_full` profile (`conf/test_full.config`): the six-run PXD024427 FragPipe LFQ tutorial dataset.
- qpx export to `<outdir>/qpx/` (on by default, `--skip_qpx_export` turns it off): bigbio's `QPX_DIANN` ([bigbio/nf-modules](https://github.com/bigbio/nf-modules) `qpx/diann`, patched) with `QPXEXPORT` (DIA-NN final quantification), `QPXFRAGPIPE` (label-free IonQuant) and `TMTANNDATA` (TMT-Integrator) write `.h5mu`, `.h5ad` and, for DIA-NN and FragPipe LFQ, the QPX Parquet dataset with [qpx](https://github.com/bigbio/qpx) 1.1.5; `--qpx_qvalue_threshold` sets the DIA-NN precursor q-value filter.

### Changed

- The DIA workflow names the mzML converted from a `.raw` file after the file (was the sample name), so DIA-NN's run names match the SDRF `comment[data file]` the qpx export checks.
- Tasks never write their inputs: FragPipe tasks copy only the staged inputs a tool misreads or writes through, and the docker tests mount inputs read-only (`docker.writableInputMounts = false` in `tests/nextflow.config`).
- nf-core template 3.5.2 -> 4.1.0.
- `manifest.nextflowVersion` `!>=25.04.0` -> `!>=25.10.4`.
- `--precursor_mass_tolerance` default 20 -> 5 ppm and `--fragment_mass_tolerance` default 20 -> 18 ppm.
  The DIA-NN workflow now passes these parameters and their units to `QUANTMSUTILS_DIANNCFG`; before, it passed `--diann_library_ms1_acc` (5) and `--diann_library_mass_acc` (18) with a fixed `ppm`, so the values DIA-NN gets by default are unchanged.
  They are Orbitrap values and equal `conf/instruments/thermo_ascend.config`; `conf/instruments/bruker_flex.config` sets 15 / 15.
- `--organism` has no default (was `Homo sapiens`): the SDRF organism comes from the `--proteome` entry, else it is guessed from the searched FASTA (UniProt `OS=` headers).
- `--fragpipe_workflow` may be omitted: `--mode fragpipe` then uses `assets/LFQ-MBR.workflow`; `--tmt_mode quant` needs one.
- `--tmt_mode` is validated against `labelcheck` and `quant`.
- FragPipe tool parameters come only from the `.workflow` file, as FragPipe reads them; the split-database slice count is `msfragger.misc.slice-db` (was `--msfragger_num_slices`).
- FragPipe groups files by the samplesheet `experiment` and `bioreplicate` columns; `condition` no longer drives grouping and only fills the TMT annotation and MSstats condition fields.
- `.raw` files reach FragPipe tools as they are; `fragpipe_convert` no longer runs ThermoRawFileParser.
- Every FragPipe module runs in the image given by `--fragpipe_container` (no default); the hardcoded `fcyucn/fragpipe:24.0` and `:24.0-academic` images are gone.
- The `test` profile is a minimal FragPipe DDA LFQ run (two human SILAC mzML subsets, `assets/LFQ-MBR.workflow`) and needs `--fragpipe_container`; it was a DIA-NN run.
- `test_tmt` reads `tests/fixtures/samplesheets/test_tmt.csv` (with `experiment` / `bioreplicate`) and adds decoys and contaminants with `--add_philosopher_decoy` / `--add_philosopher_contaminants`.
- The DIA experimental design is built from the prepared spectra files and lists every sample with a numeric sample index and its `bioreplicate`.
- `outputDir = params.outdir` and `workflow.output.mode = params.publish_dir_mode` are set in `nextflow.config` (template 4.1.0); outputs are still published with `publishDir`.
- `env.LANG = 'C.UTF-8'`: FragPipe's Java tools run with UTF-8 as desktop FragPipe does (with the image's default ASCII encoding the matched-fragment ion order differed).
- The pipeline-wide `errorStrategy = 'retry'`, `maxRetries` and `maxErrors` are removed from `nextflow.config`; tasks retry only on the `conf/base.config` exit codes (130-145, 104, 175-177; 176 and 177 are new), and `DECOMPRESS` / `TDF2MZML` retry always.
- `process.resourceLimits` moved to `conf/base.config` (96 CPUs, 760 GB, 72 h), loaded before the profiles so test profiles can lower them.
- Local modules are renamed to nf-core 4.1.0 module names (alphanumeric directory, process = upper-cased directory):

  | Old                                                                 | New                                                           |
  | ------------------------------------------------------------------- | ------------------------------------------------------------- |
  | `fragpipe_headless` / `FRAGPIPE_HEADLESS`                           | `fragpipeheadless` / `FRAGPIPEHEADLESS`                       |
  | `generate_sdrf_from_samplesheet` / `GENERATE_SDRF_FROM_SAMPLESHEET` | `generatesdrffromsamplesheet` / `GENERATESDRFFROMSAMPLESHEET` |
  | `merge_split_search` / `MERGE_SPLIT_SEARCH`                         | `mergesplitsearch` / `MERGESPLITSEARCH`                       |
  | `msfragger_calibrate` / `MSFRAGGER_CALIBRATE`                       | `msfraggercalibrate` / `MSFRAGGERCALIBRATE`                   |
  | `msstats_lfq` / `MSSTATS_LFQ`                                       | `msstatslfq` / `MSSTATSLFQ`                                   |
  | `mzml_indexing` / `MZML_INDEXING`                                   | `mzmlindexing` / `MZMLINDEXING`                               |
  | `openms_peak_picker` / `OPENMS_PEAK_PICKER`                         | `openmspeakpicker` / `OPENMSPEAKPICKER`                       |
  | `parse_fragpipe_workflow` / `PARSE_FRAGPIPE_WORKFLOW`               | `parsefragpipeworkflow` / `PARSEFRAGPIPEWORKFLOW`             |
  | `philosopher/database` / `PHILOSOPHER_DATABASE`                     | `philosopherdatabase` / `PHILOSOPHERDATABASE`                 |
  | `philosopher/filter` / `PHILOSOPHER_FILTER`                         | `philosopherfilter` / `PHILOSOPHERFILTER`                     |
  | `philosopher/labelquant` / `PHILOSOPHER_LABELQUANT`                 | `philosopherlabelquant` / `PHILOSOPHERLABELQUANT`             |
  | `philosopher/report` / `PHILOSOPHER_REPORT`                         | `philosopherreport` / `PHILOSOPHERREPORT`                     |
  | `split_fasta` / `SPLIT_FASTA`                                       | `splitfasta` / `SPLITFASTA`                                   |
  | `tmtlabelcheck` / `TMT_LABELCHECK_ANALYZE`                          | `tmtlabelcheckanalyze` / `TMTLABELCHECKANALYZE`               |

- `bin/generate_sdrf.py` moved out of the `generate_sdrf_from_samplesheet` module's `resources/usr/bin/`.
- Patched nf-core components (`modules.json` `patch`):
  - `diann`: `predict_speclib` matches `*.predicted.speclib` (DIA-NN names it after the FASTA, not `--out-lib`); the stub writes the command line and a settings line in the log.
  - `thermorawfileparser`: the stub reports a fixed version instead of running the tool (new patch).
  - `dia_proteomics_analysis`: no `null` channel items (empty optional inputs are `[[]]`, presence is tested with `count()`); emits `pg_matrix` and `diann_log`.
  - `utils_nextflow_pipeline`: `groovy.json.JsonGenerator$Options` (new patch).
  - `utils_nfcore_pipeline`: patch removed (unpatched nf-core version).
- `.nf-core.yml` pins `diann` (DIA-NN 1.8.1) and `utils_nfschema_plugin` (nf-schema 2.5.1) for `nf-core modules update`.
- `validation.defaultIgnoreParams` is `["iproteomes"]` (was `["genomes", "databases"]`).
- `CONTRIBUTING.md` moved from `.github/` to `docs/`.

### Fixed

- FragPipe mode no longer publishes Philosopher's `.meta` workspace (`ionquant/`, `opair/`, `mbg/`, `ptmshepherd/groups/`, `freequant/`, `labelquant/`): as FragPipe 24.0 does (FragpipeRun.java:2544-2583), the new `PHILOSOPHERWORKSPACECLEAN` runs `philosopher workspace --clean --nocheck` in each group dir after Filter's report, and the FreeQuant and LabelQuant workspaces are no longer published.
- `--diann_skip_preliminary_analysis` reads `--diann_empirical_assembly_lib` (it read `--diann_speclib`) and stops when the library or `--diann_empirical_assembly_log` is missing.
- `dia_proteomics_analysis`: a given empirical library now skips preliminary analysis and library assembly; the `join` on the FASTA-keyed channel never matched the experiment-keyed library, so they always ran.
- DIA experimental design: samples whose `condition` was empty or equal to the sample name were left out.

### Removed

- `DIANNR` (DIA-NN R matrices); `IQ` gives the MaxLFQ tables.
- `SKYLINE` (FragPipe on Linux does not run Skyline).
- `MSFRAGGER_INDEX`: MSFragger builds its pepindex inside every search task, as FragPipe does.
- Unused modules and subworkflows: `FRAGPIPE` (all-in-one), `PHILOSOPHER_PIPELINE`, `MZML_STATISTICS` (the nf-core `QUANTMSUTILS_MZMLSTATISTICS` runs), `SAMPLESHEET_CHECK`, `SDRF_PARSING`, `PREPROCESS_EXPDESIGN`, `FASTQC`, `INPUT_CHECK`, `CREATE_INPUT_CHANNEL`; SDRF input was never wired, the CSV samplesheet stays the only input.
- `lib/WorkflowUtils.groovy`; its functions are Nextflow functions in `subworkflows/local/utils_nfcore_msproteomics_pipeline/`.
- `conf/reference_proteomes.config` and `conf/reference_proteomes_ignored.config` (UniProt download), replaced by iProteomes.
- `conf/igenomes.config` and `conf/igenomes_ignored.config` (`skip_features: igenomes`).
- `conf/tests/test.config`.
- Slack / Microsoft Teams notifications (`--hook_url`, `assets/adaptivecard.json`, `assets/slackreport.json`).
- `.github/workflows/linting_comment.yml`.

### Known issues

- `params.pipelines_testdata_base_path` points at the `an-altosian/test-datasets` fork (`msproteomics` branch) until nf-core/test-datasets#1953 merges.
  After the merge it goes back to `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/msproteomics`, together with the spectra URLs in `tests/fixtures/samplesheets/test_*.csv`.

### Parameters

| Old parameter                           | New parameter                    |
| --------------------------------------- | -------------------------------- |
| `--proteomes_base`                      | `--iproteomes_base`              |
| `--proteomes_ignore`                    | `--iproteomes_ignore`            |
|                                         | `--proteome`                     |
|                                         | `--save_reference`               |
|                                         | `--add_philosopher_decoy`        |
|                                         | `--add_philosopher_contaminants` |
|                                         | `--fragpipe_threads`             |
|                                         | `--metaproteomics_names_dmp`     |
|                                         | `--metaproteomics_nodes_dmp`     |
| `--acquisition_method`                  |                                  |
| `--add_decoys`                          |                                  |
| `--config_file`                         |                                  |
| `--decoy_method`                        |                                  |
| `--decoy_string`                        |                                  |
| `--decoy_string_position`               |                                  |
| `--decoy_tag`                           |                                  |
| `--enable_msstats`                      |                                  |
| `--fragment_method`                     |                                  |
| `--genome`                              |                                  |
| `--hook_url`                            |                                  |
| `--id_only`                             |                                  |
| `--is_phospho_data`                     |                                  |
| `--labelling_type`                      |                                  |
| `--local_input_type`                    |                                  |
| `--msfragger_num_slices`                |                                  |
| `--mzml_features`                       |                                  |
| `--peptide_fdr`                         |                                  |
| `--philosopher_database`                |                                  |
| `--protein_fdr`                         |                                  |
| `--psm_fdr`                             |                                  |
| `--shuffle_max_attempts`                |                                  |
| `--shuffle_sequence_identity_threshold` |                                  |
| `--skip_experimental_design_validation` |                                  |
| `--skip_factor_validation`              |                                  |
| `--skip_ms_validation`                  |                                  |
| `--test_mode`                           |                                  |
| `--tmt_annotation`                      |                                  |
| `--use_ols_cache_only`                  |                                  |
| `--validate_ontologies`                 |                                  |

Changed defaults:

| Parameter                        | Old default                                                | New default                                                                           |
| -------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `--precursor_mass_tolerance`     | 20                                                         | 5                                                                                     |
| `--fragment_mass_tolerance`      | 20                                                         | 18                                                                                    |
| `--organism`                     | `Homo sapiens`                                             | none                                                                                  |
| `--pipelines_testdata_base_path` | `https://raw.githubusercontent.com/nf-core/test-datasets/` | `https://raw.githubusercontent.com/an-altosian/test-datasets/refs/heads/msproteomics` |

### Software dependencies

| Dependency | Old version                                 | New version                   |
| ---------- | ------------------------------------------- | ----------------------------- |
| MultiQC    | 1.33                                        | 1.35                          |
| OpenMS     | `bigbio/openms-tools-thirdparty` 2025.04.14 | 3.5.0                         |
| IQ (R)     | `dongzehe/iq:0.0.1`                         | 2.0.1                         |
| FragPipe   | 24.0 (fixed images)                         | 24.0 (`--fragpipe_container`) |
| DIA-NN     | 1.8.1                                       | 1.8.1 (pinned)                |
| nf-schema  | 2.5.1                                       | 2.5.1                         |

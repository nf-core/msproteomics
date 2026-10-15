---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/msproteomics`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/msproteomics.

This page describes the recommended nf-core way to contribute to both nf-core/msproteomics and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/msproteomics.

> [!NOTE]
> If you need help using or modifying nf-core/msproteomics, ask on the nf-core Slack [#msproteomics](https://nfcore.slack.com/channels/msproteomics) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/msproteomics/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/msproteomics repository](https://github.com/nf-core/msproteomics) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focused as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/msproteomics documentation](https://nf-co.re/msproteomics/usage) or ask on the nf-core [#msproteomics](https://nfcore.slack.com/channels/msproteomics) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/msproteomics without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/msproteomics/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag stub --profile +docker --verbose   # what CI runs (make test)
make test-python                                     # Python unit tests (tests/python/requirements.txt)
```

Tests tagged `integration` need real data and, for FragPipe modules, the licensed FragPipe image passed at runtime; golden tests (`tests/golden*.nf.test`, tag `golden`) are only collected when `GOLDEN_ROOT` is set (see [Golden FragPipe tests](#golden-fragpipe-tests)).

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files.
Every test asserts a snapshot of its deterministic output, except where the test compares with something other than a recorded snapshot: golden tests compare with FragPipe's own outputs, and failure tests assert the error message.

Some modules have stub tests only, because real inputs are private data or need the licensed image:
`opair`, `pairscans`, `philosopherabacus`, `philosopherfilter`, `philosopheriprophet`, `philosopherlabelquant`, `philosopherproteinprophet`, `philosopherreport`, `ptmprophet`, `ptmshepherd`, `tmtintegrator` and `writesubmzml`.
Their real-data coverage is the golden tests (`modules/local/<module>/tests/golden*.nf.test`, compared with FragPipe 24.0 runs on the same data) and the Tower parity runs of the whole pipeline.
`pmultiqc` runs for real in the DIA pipeline test (`tests/default.nf.test`, tag `dia_integration`).
`openmspeakpicker` (`--openms_peakpicking`) and `tdf2mzml` (`--convert_dotd`) are off by default and have no real-data test yet.
Update the snapshots with the following command:

```bash
nf-test test --tag stub --profile +docker --verbose --update-snapshot
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/msproteomics code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

- FragPipe tool parameters live in `.workflow` files (`assets/*.workflow`), not in `nextflow.config`.
- Modules only run tools: argv strings and config files are built by `parse_fragpipe_workflow.py` (`modules/local/parsefragpipeworkflow/`) or the Groovy helpers in `subworkflows/local/utils_nfcore_msproteomics_pipeline/`, each with a unit test.
  One approved exception: the FragPipe jar lookup (`fragpipeJarLines`, `subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a.nf`) runs `find`, `awk`, `sort` and `unzip` in the task.
  The licensed FragPipe image is a runtime-only parameter, so the MSFragger, IonQuant and diaTracer jar names (`MSFragger-4.4.1.jar`, or with an edition word before the version) and the batmass-io and MSBooster jars the image's FragPipe build names are only known inside the task.
  The lookup copies FragPipe's own: `TabConfig.getJarPath` (TabConfig.java:176-178, :188-209, the highest `<tool>-<version>.jar` under the tools folder), MSFragger's `ext/` libraries next to its jar (FragpipeRun.java:1348-1356), and the jar names compiled into `ToolingUtils.BATMASS_IO_JAR` (ToolingUtils.java:62) and `CmdMSBooster.JAR_MSBOOSTER_NAME` (CmdMSBooster.java:52).
- Local modules import those helpers with `include { ... } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/...'`.
  nf-core modules are self-contained; this pipeline keeps the helpers (about 4,000 lines of FragPipe argv and file builders, shared by many modules) in one tested place instead of copying them into each module.
- `ext.args` holds raw command-line arguments only (TMT-Integrator, which has no CLI arguments, is the exception).
- The FragPipe container image is a runtime parameter (`--fragpipe_container`); do not name an image in the repository.
- Every new module or subworkflow needs a stub nf-test with a snapshot; run `make lint` and `make test` before opening a pull request.

## Golden FragPipe tests

Golden tests check that a module reproduces FragPipe 24.0 exactly: the module gets the exact input files FragPipe used, and its command line and every output file are compared with FragPipe's.
They need licensed-environment data (a captured FragPipe run and the licensed FragPipe image), so they run outside public CI with the nf-test tag `golden`.

### Capture a FragPipe run

Run FragPipe 24.0 headless, keeping every intermediate file, then capture it:

```bash
python3 tests/tools/golden_capture.py \
    --log <run>/work/log_<timestamp>.txt \
    --dataset <dataset name> \
    --workflow <workflow name> \
    --root /path/to/golden
```

If the run dir holds empty placeholder spectra files (for example Docker bind-mount targets), pass `--lcms-source <dir with the real files>`; an empty spectra file without it is an error.
`--root` defaults to `$GOLDEN_ROOT` and is required when that is unset; an unwritable root is an error.
The capture writes `<root>/<dataset>/<workflow>/` with:

| Path                                       | Content                                                                                                                                                                                                                                                                   |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `manifest.json`                            | every executed FragPipe step: tool, normalized argv, work dir, exit code, parallel batch and concurrent steps, inputs (with `golden_state`), outputs (with final state: present, moved, deleted, overwritten), flags for anything inferred; sha256 of every captured file |
| `fragpipe_log.redacted.txt`                | the FragPipe log with license arguments and license lines removed                                                                                                                                                                                                         |
| `files/work/`, `files/input/`, `files/db/` | FragPipe work dir, spectra dir (symlinks resolved), database dir                                                                                                                                                                                                          |

The FragPipe log is never copied verbatim, license-named files are refused, and the capture aborts if a copied file references a license file or `--license` argument.
In the redacted log and the manifest, the value of any license flag (`--license`, `-licence`, `--License-key`, ..., with a space or `=`) becomes `{LICENSE}` in place and any `license.dat` path becomes `{LICENSE_PATH}`; the argument count is unchanged, command lines are never dropped, and only informational lines that mention a license are dropped.
Logs already scrubbed this way are accepted.

FragPipe runs the per-group steps of a parallel group concurrently, so their log output interleaves.
The manifest records each step's `parallel_batch` and `concurrent_with`, and flags steps whose `Process '<name>' finished` marker cannot be told apart from a concurrent step of the same name.

Some steps rewrite their inputs in place (IonQuant, PTM-Shepherd in glyco mode, MBG, `Rewrite pepxml`, `AppendFile`), so a captured run only holds the rewritten version.
Each input carries a `golden_state` (`available`, `fragpipe_install`, `rewritten_by_step:N`, `moved_by_step:N`, `deleted_by_step:N`), and `golden_inputs_available` is false when any input is gone.
Golden tests for such steps take their inputs from a run that stops just before the step, for example a run that ends after Philosopher Report.
`tests/tools/fragpipe_log_to_manifest.py <log> -o manifest.json` builds the manifest alone.

### Run the golden tests

```bash
export FRAGPIPE_CONTAINER=<licensed FragPipe 24.0 image>   # the FragPipe tools of the golden tests run in it
export GOLDEN_ROOT=/path/to/golden                         # captured FragPipe runs; nf-test.config collects the golden tests only when it is set
export OPAIR_DET_DLL=/path/to/PTMLocalization.dll          # deterministic O-Pair 1.1.10 build (O-Pair golden test)
NXF_OFFLINE=true nf-test test --tag golden
```

Stock O-Pair 1.1.10 breaks ties between equally scored glycan placements in a random order, so its `psm.tsv` differs between identical runs.
The O-Pair golden test compares against `PXD017646/glyco-O-Pair_perfile_det`, a FragPipe run with a `PTMLocalization.dll` that builds its glycan boxes in enumeration order (the order of stock O-Pair with `-t 1`).
It mounts `$OPAIR_DET_DLL` over the image's `tools/opair/PTMLocalization.dll`, checks its md5 (`01fca5e12c419a45880c5b610bcd0175`) and requires `psm.tsv` to match exactly.

Each golden test has its own `golden.config` next to it: it runs the tools in `$FRAGPIPE_CONTAINER` with docker and pins the tested processes to the threads and memory of the golden FragPipe run.

Each test (`modules/**/tests/golden.nf.test`) selects its steps from the manifest, stages the golden inputs, runs the module, and calls `tests/tools/golden_step_check.py`, which fails on any of:

- argv differences after path normalization (step work dir, FragPipe install and external inputs become placeholders; argument order and flags are compared as is);
- output files that differ under the comparison rules below, or that the module did not write;
- module files with no FragPipe counterpart, except module logs passed with `--ignore-actual`.

Outputs FragPipe deleted or overwrote later (for example Percolator's `*_psms.tsv`) have no golden copy and are reported as `UNVERIFIABLE`; a step with no output that could be compared fails (`NO_VERIFIED_OUTPUT`).
The module's command line is split with POSIX shell rules, so quoted and empty arguments are compared like FragPipe's.

### Comparison rules

`tests/tools/golden_compare.py` applies exactly the rules measured by running FragPipe 24.0 twice on the same input; every other byte must match.
The Philosopher tolerances are the measured envelope of its nondeterministic FDR tie handling (`tests/tools/philosopher_tables.py`, determinism doc section "Philosopher filter/report nondeterminism").
There are two rule profiles, selected with `--profile` (default `dda`); `tests/tools/golden_step_check.py` picks `dia` when the golden run's spectra are DIA.

Profile `dda` (TMT label check, DDA LFQ, glyco DDA):

| File                                                                                                     | Normalization                                                                                                                      |
| -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `log_*.txt`, anything under `?/`                                                                         | skipped                                                                                                                            |
| `*.pepindex`                                                                                             | first 8 bytes dropped, run paths rewritten                                                                                         |
| `*.prot.xml`                                                                                             | `<tmpdir>/<uuid>/`, run paths, ISO `date=`/`time=`                                                                                 |
| `*.pepXML`, `*.pep.xml`                                                                                  | run paths, ISO `date=`/`time=`                                                                                                     |
| `filter.log`                                                                                             | run paths, `time=` stripped, `Parsing` lines as a set, PSM/peptide/ion counts ±2, decoy ±1, threshold ±2e-4 (Philosopher envelope) |
| `fragger*.params`, `*.workflow`, `*.job`, `*.fp-manifest`, `filelist_*.txt`, `experiment_annotation.tsv` | run paths                                                                                                                          |
| `psm.tsv`                                                                                                | Philosopher envelope: `Qvalue` within 1.01/max(T_above,1), at most 1 boundary PSM between the two runs' `filter.log` thresholds    |
| `peptide.tsv`, `ion.tsv`, `protein.tsv`                                                                  | spectral counts may move by the boundary PSMs mapped to the entry                                                                  |
| everything else                                                                                          | none (raw bytes)                                                                                                                   |

Profile `dia` adds the rules measured on the DIA runs:

| File                                             | Normalization                                                                                                                                                                           |
| ------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MSBooster/spectraRT.predicted.bin`              | skipped                                                                                                                                                                                 |
| `*.pdf` (EasyPQP alignment plots)                | `/CreationDate (D:...)`                                                                                                                                                                 |
| DIA-NN `report.log.txt`                          | run paths, `Current date and time:` text, `[m:ss] ` prefixes                                                                                                                            |
| `MSBooster/spectraRT*.tsv`                       | rows sorted (header first)                                                                                                                                                              |
| `peptide.tsv`, `ion.tsv`                         | full Philosopher envelope: `Qvalue` rule; `ion.tsv` keyed on stripped sequence + `Charge` + `M/Z`, `Expectation`/`Prev AA`/`Next AA` dropped; `peptide.tsv` without `Prev AA`/`Next AA` |
| `msbooster_params.txt`, `easypqp_files/cmds.txt` | run paths                                                                                                                                                                               |

Compare two run directories, for example two FragPipe runs:

```bash
python3 tests/tools/golden_compare.py dirs runA/work runB/work --profile dda \
    --expected-prefix runA --actual-prefix runB \
    --expected-tmp runA/jtmp --actual-tmp runB/jtmp
```

Missing and extra files are failures, and so is a comparison in which no file was compared.

### Coverage

```bash
python3 tests/tools/golden_coverage.py <golden>/<dataset>/<workflow>/manifest.json
```

lists every FragPipe step of the run and the golden test that covers it.
Unit tests for the scripts: `make test-python` (dependencies in `tests/python/requirements.txt`).

## Comparing against FragPipe

`tests/tools/compare_fragpipe_results.py` measures how closely an msproteomics run reproduces a FragPipe GUI/headless run on the same data.
It needs only Python 3 (standard library).

```bash
python3 tests/tools/compare_fragpipe_results.py \
    --test results/ \
    --reference /path/to/fragpipe_output/ \
    --outdir comparison/
```

Both directories are searched recursively, so the msproteomics layout (`philosopher/<exp>/`, `ionquant/`, `tmtintegrator/`) and the FragPipe layout (`<exp>/`, `combined_*.tsv`, `tmt-report/`) can be compared directly.
Per-experiment `psm.tsv`, `peptide.tsv`, `ion.tsv` and `protein.tsv` are matched by experiment directory name.
Use `--experiment-map TEST=REF ...` when experiment names differ, and `--exclude REGEX` when a file is found twice.
Combined tables (`combined_protein.tsv`, `combined_peptide.tsv`, `combined_modified_peptide.tsv`, `combined_ion.tsv`) and TMT-Integrator tables (`abundance_*.tsv`, `ratio_*.tsv`) are matched by file name.
A table present on only one side, or a missing required column, is an error (exit code 2) that names the file.

The tool reports, per experiment and for the combined tables:

| Metric         | Description                                                                                                                 |
| -------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Counts         | PSMs, peptides, modified peptides, ions and proteins (decoys excluded)                                                      |
| Overlap        | Jaccard, only-test and only-reference counts for peptide sequences, modified peptides, ions and protein IDs                 |
| Shared PSMs    | Percentage of shared spectra with the same peptide and modified peptide, and Probability/Hyperscore/Expectation correlation |
| Quantification | Pearson and Spearman of log2 intensities (raw and MaxLFQ), median log2 ratio and percentage within 2-fold                   |
| TMT            | Per-channel correlation of TMT-Integrator log2 values                                                                       |
| Parameters     | Differing keys between `fragger.params` and `fragpipe.workflow` (informational, not gated)                                  |

Modified peptides are keyed on `Assigned Modifications` (position, residue and mass rounded to `--mod-decimals`), because `Modified Peptide` notation differs between tables.
The verdict is PASS (exit code 0) only when every gated metric meets `--max-count-deviation` (default 0.01), `--min-jaccard` (default 0.98) and `--min-spearman` (default 0.99); otherwise it is FAIL (exit code 1).
Correlations with fewer than `--min-quant-n` shared values are not computed and are listed as skipped.
Results are written to `comparison.md` and `comparison.json` in `--outdir`.

Unit tests run with `make test-python`.

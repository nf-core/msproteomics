# PLANTOTOOLCONFIGS

Tool configs for the FragPipe subworkflows from FragPipe's configured run.

## Description

PLANTOTOOLCONFIGS runs `bin/plan_to_tool_configs.py` on FRAGPIPEPLAN's work dir (`--write-plan`).
The steps come from the run log (`fragpipe_log_to_manifest.py`: `parse_dry_run`, `decode_command`, `step_kind`), each step's config files are the work-dir files its argv names, and paths go under the run's roots (`{WORKDIR}`, `{LCMS_DIR}`, `{DB_DIR}`, `{FRAGPIPE}`, the parent of FragPipe's tools folder, `{TMP}`); the licence becomes `{LICENSE}`.
When several argvs print as a command line, the one whose differing arguments are the run's values (the saved `fragpipe.workflow` values, the LC-MS, database and work-dir paths, optionally in FragPipe's literal quotes) is taken; a line with no such reading or two-plus, and a masked `--api-key`, stop the script.
It writes the outputs the subworkflows read: `tool_configs.json` (run flag, args and fields per tool), the MSFragger params of each data type the run has, and the MSBooster params.
Each tool's `args` is the step's argv minus what the msproteomics module builds itself (binary, JVM options, threads, licence, inputs and outputs); params files are FragPipe's own bytes.
The script also stops runs FragPipe plans but msproteomics cannot run: MSFragger off or digest-only, FreeQuant without the Philosopher filter, psm.tsv readers without the Philosopher report, transfer learning, Write sub mzML with a split database, IonQuant values with whitespace.
A plan step with no extractor, per-file steps with different options, or a feature FragPipe cannot plan before the run (two-pass `_sub.mzML`) stops the script with an error.

## Inputs

| Channel | Type     | Description                        |
| ------- | -------- | ---------------------------------- |
| meta    | val(map) | Plan id map                        |
| workdir | path     | FragPipe's work dir (FRAGPIPEPLAN) |

## Outputs

| Channel                   | Type | Description                                                         |
| ------------------------- | ---- | ------------------------------------------------------------------- |
| tool_configs_json         | path | `tool_configs.json`                                                 |
| msfragger_config          | path | `configs/fragger.params` or `configs/fragger_dda.params` (optional) |
| msfragger_dia_config      | path | `configs/fragger_dia.params` (optional)                             |
| msfragger_gpfdia_config   | path | `configs/fragger_gpfdia.params` (optional)                          |
| msfragger_dda_plus_config | path | `configs/fragger_dda_plus.params` (optional)                        |
| msbooster_config          | path | `configs/msbooster.config` (empty when MSBooster does not run)      |
| versions_python           | val  | Python version (`versions` topic)                                   |

## Container

`quay.io/biocontainers/python:3.12`.

## Tests

`tests/python/test_plan_to_tool_configs.py` (`PlannerParityTest`) compares the tool configs with the former parser's tool configs, frozen next to each golden fixture (`tests/fixtures/fragpipe_plans/<golden>.parser_tool_configs.json`), on the golden work dirs in `tests/fixtures/fragpipe_plans/<golden>/`.

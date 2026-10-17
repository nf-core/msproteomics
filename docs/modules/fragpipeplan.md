# FRAGPIPEPLAN

FragPipe's configured run for a workflow and manifest, without running any tool.

## Description

FRAGPIPEPLAN runs FragPipe headless with `--write-plan --msfragger-per-file`.
FragPipe configures every tool as a real run would and writes their config files, then stops before starting any process.
Its work dir `fp_work/` holds those config files, FragPipe's run files (`fragpipe.workflow`, `fragpipe-files.fp-manifest`, `fragpipe.job`) and the run log `log_<time>.txt`, which lists every command in order with its work dir and command line.

The task stages the LC-MS files as links (FragPipe only checks that they exist) and the FASTA (FragPipe counts its entries and decoys).
It writes the manifest and the run's `.workflow` with the FRAGPIPE_HEADLESS helpers (`headlessLayout`, `headlessStagingLines`, `headlessWorkflowText`).
FragPipe creates its work dir `fp_work/`.
The portable copies (`bin/fragpipe_portable.py`) have every staged input at its original path, and the log has its licence material redacted.
MSFragger always plans its per-file search (`--msfragger-per-file`), as msproteomics searches each file in its own MSFRAGGER task.

## Inputs

| Channel         | Type        | Description                                                                       |
| --------------- | ----------- | --------------------------------------------------------------------------------- |
| meta            | val(map)    | Plan id map                                                                       |
| raw_files       | path        | LC-MS files, staged into `raw_files/` as links                                    |
| database        | path        | Searched FASTA, with decoys                                                       |
| workflow_text   | val(string) | `.workflow` text                                                                  |
| manifest_rows   | val(list)   | `[file name, experiment, bioreplicate, data type, FragPipe group]` per LC-MS file |
| annotation_rows | val(list)   | `[group, channel, sample name]` per TMT channel, `[]` for label-free              |

## Outputs

| Channel           | Type | Description                                                          |
| ----------------- | ---- | -------------------------------------------------------------------- |
| workdir           | path | `fp_work/`: run log, tool config files, FragPipe's run files         |
| manifest          | path | The manifest FragPipe read                                           |
| portable_log      | path | The run log, licence redacted, inputs at their original paths        |
| portable_workflow | path | FragPipe's saved `fragpipe.workflow`, inputs at their original paths |
| portable_manifest | path | FragPipe's saved manifest, inputs at their original paths            |
| versions          | path | Software versions (`versions.yml`)                                   |

## Parameters

| Parameter  | Default                                   | Description                                                                  |
| ---------- | ----------------------------------------- | ---------------------------------------------------------------------------- |
| `ext.args` | `--threads <--fragpipe_threads>` when set | FragPipe CLI arguments; without `--threads` FragPipe uses `workflow.threads` |
| (memory)   | 90% of task.memory                        | FragPipe's `-Xmx` and `--ram`                                                |

## Container

A FragPipe build with the plan options proposed in [Nesvilab/FragPipe#2907](https://github.com/Nesvilab/FragPipe/pull/2907) (`--write-plan` writes the tool config files and the commands without running them; `--msfragger-per-file`), passed with `--fragpipe_container`.
No FragPipe release has it yet; with any other FragPipe the task stops with an error naming the missing option.

## Usage

Run by FRAGPIPE_TOOL_CONFIGS (not with `--fragpipe_plan`, which gives the work dir instead); PLANTOTOOLCONFIGS reads its work dir.

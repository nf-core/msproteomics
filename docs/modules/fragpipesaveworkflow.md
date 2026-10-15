# FRAGPIPESAVEWORKFLOW

Write FragPipe's saved UI state (`fragpipe.workflow`) for a run.

## Description

FRAGPIPESAVEWORKFLOW writes the `fragpipe.workflow` file that FragPipe writes into its work directory before any tool runs (`FragpipeRun.java:371-375`).
`assets/fragpipe_ui_state/FragpipeUiState` drives FragPipe's own headless load and save code, so the file holds every UI key with its default, the run values and the tool-version header.
It needs the FragPipe image.

## Inputs

| Channel                                                             | Type  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ------------------------------------------------------------------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, workflow_file, manifest, ui_state, database, workdir, spectra | tuple | `workflow_file`: the run's `.workflow` file staged under `input/` (so a file named `fragpipe.workflow` is not overwritten through its staged link); `manifest`: `fragpipe-files.fp-manifest` naming the staged spectra; `ui_state`: directory with the compiled `FragpipeUiState` class; `database`: the run database, written as `database.db-path`; `workdir`: absolute path written as `workdir`; `spectra`: the LC-MS files named in the manifest (FragPipe's headless manifest load stops when a listed file does not exist) |

## Outputs

| Channel           | Type                    | Description         |
| ----------------- | ----------------------- | ------------------- |
| workflow          | tuple(val, path)        | `fragpipe.workflow` |
| versions_fragpipe | tuple (topic: versions) | FragPipe version    |

## Parameters

| Parameter | Default | Description                 |
| --------- | ------- | --------------------------- |
| ext.args  | ''      | Additional CLI arguments    |
| ext.when  | (none)  | Skip the process when false |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in `workflows/fragpipe.nf` and `workflows/tmt_labelcheck.nf`; its output feeds FRAGPIPESDRF and FRAGPIPESECONDPASS.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

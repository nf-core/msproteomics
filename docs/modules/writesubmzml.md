# WRITESUBMZML

FragPipe "Write sub mzML": the scans no confident PSM explains.

## Description

WRITESUBMZML mirrors `CmdWriteSubMzml` / `util.WriteSubMzml` (FragPipe 24.0, `FragpipeRun.java:2388-2404`, `CmdWriteSubMzml.java:97-130`).
It runs per LC-MS file after the final Philosopher report when `tab-run.write_sub_mzml=true`.
WriteSubMzml reads `<stem>.mzBIN_all` next to the LC-MS file, drops the MS2 scans whose PSM probability in `psm.tsv` is above the threshold (and the MS3 scans of their parents), writes the other scans with MSFragger's calibrated peaks into `<stem>_sub.mzML` and deletes the `.mzBIN_all` on exit.
The argv comes from `writeSubMzmlArgv` (`fragpipe_utils_a.nf`).

## Inputs

| Channel                                    | Type  | Description                                                                                                                                                                                                                                |
| ------------------------------------------ | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| meta, lcms, mzbin_all, psm, prob_threshold | tuple | `meta.group` is the file's FragPipe group; `lcms`: the LC-MS file as given; `mzbin_all`: MSFragger's `<stem>.mzBIN_all` (`write_mzbin_all = 1`); `psm`: `psm.tsv` of the file's group; `prob_threshold`: `tab-run.sub_mzml_prob_threshold` |

## Outputs

| Channel           | Type                    | Description                                                        |
| ----------------- | ----------------------- | ------------------------------------------------------------------ |
| sub_mzml          | tuple(val, path)        | `<stem>_sub.mzML`                                                  |
| log               | tuple(val, path)        | `write_sub_mzml.log`, the command line followed by the tool output |
| versions_fragpipe | tuple (topic: versions) | FragPipe version                                                   |

## Parameters

| Parameter | Default | Description                 |
| --------- | ------- | --------------------------- |
| ext.args  | ''      | Additional CLI arguments    |
| ext.when  | (none)  | Skip the process when false |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in `workflows/fragpipe.nf` and the `fragpipe`, `fragpipe_search` and `tmt_labelcheck` subworkflows; its output feeds a second run (see FRAGPIPESECONDPASS).

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

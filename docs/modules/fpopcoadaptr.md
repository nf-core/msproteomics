# FPOPCOADAPTR

FPOP coADAPTr converter, run once over the whole FragPipe work dir.

## Description

FPOPCOADAPTR mirrors FragPipe 24.0 `CmdFPOPcoadaptr.java`.
It edits the final `psm.tsv` files (or the DIA-NN `report.tsv` for the `diann` layout, where `FPOP_DiannReport_Editor` runs instead of `FPOP_PSMwriter`, `CmdFPOPcoadaptr.java:429-433`) and writes the result under the output prefix.

## Inputs

| Channel                         | Type                       | Description                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, inputs, layout, fpop_mods | tuple(val, path, val, val) | `inputs`: what FragPipe's `Files.walk` finds in the work dir, staged under `input/<n>/`; `layout`: `tree` (IonQuant results directories, copied whole), `groups` (Philosopher per-group directories, only `psm.tsv` is copied; `workdir-root` is the blank group) or `diann` (DIA-NN `dia-quant-output` directories, only `report.tsv`); `fpop_mods`: comma-joined `%.4f` FPOP masses (`FpopCoadaptrPanel.getFpopMods()`) |

## Outputs

| Channel           | Type                    | Description                                                                                                                 |
| ----------------- | ----------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| tables            | tuple(val, path)        | All files of the edited result tree (`<prefix>/**`), emitted as files so they can be published at the replaced task's paths |
| psms              | tuple(val, path)        | `<prefix>/**/psm.tsv` (optional)                                                                                            |
| diann_report      | tuple(val, path)        | `<prefix>/**/report.tsv` (optional)                                                                                         |
| log               | tuple(val, path)        | `fpop_coadaptr.log`                                                                                                         |
| versions_fragpipe | tuple (topic: versions) | FragPipe version                                                                                                            |

## Parameters

| Parameter  | Default | Description                                 |
| ---------- | ------- | ------------------------------------------- |
| ext.args   | ''      | Additional CLI arguments                    |
| ext.prefix | meta.id | Override the output directory / file prefix |
| ext.when   | (none)  | Skip the process when false                 |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`).

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

# REPRINTTOSAINTINPUT

SAINTexpress inputs from a REPRINT table.

## Description

REPRINTTOSAINTINPUT ports `CmdSaintExpress.reprintTsvToSaintInput` and its input check (FragPipe 24.0, `CmdSaintExpress.java:56-81`, `:98-167`).
It writes `inter.dat`, `bait.dat` and `prey.dat`, or `skipped_mode.txt` with FragPipe's error when the table is missing (the mode is then skipped).

## Inputs

| Channel                 | Type                  | Description                                                                                                                   |
| ----------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| meta, reprint_tsv, mode | tuple(val, path, val) | `reprint.<mode>.tsv` of the work dir (Abacus `--reprint` or IonQuant), or `[]` when it does not exist; `mode`: `spc` or `int` |

## Outputs

| Channel         | Type                    | Description                                                         |
| --------------- | ----------------------- | ------------------------------------------------------------------- |
| saint_input     | tuple(val, path, val)   | `saintexpress-<mode>/{inter,bait,prey}.dat` and the mode (optional) |
| skipped         | tuple(val, path)        | `saintexpress-<mode>/skipped_mode.txt` (optional)                   |
| versions_python | tuple (topic: versions) | Python version                                                      |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`) before SAINTEXPRESS.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

# SORTFASTA

Sort the proteins of a FASTA file for the calibration search of a split-database run.

## Description

SORTFASTA runs `sort_fasta.py`, the port of `sample_fasta()` in FragPipe 24.0 `msfragger_pep_split.py` (lines 614-620, 632-633).
The proteins are sorted bytewise.

## Inputs

| Channel | Type | Description      |
| ------- | ---- | ---------------- |
| fasta   | path | Protein database |

## Outputs

| Channel         | Type                    | Description           |
| --------------- | ----------------------- | --------------------- |
| fasta           | path                    | `sorted_<fasta name>` |
| versions_python | tuple (topic: versions) | Python version        |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.args                              | ''      | Additional CLI arguments              |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`); the sorted FASTA is the database of MSFRAGGERCALIBRATE.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

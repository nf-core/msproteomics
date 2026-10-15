# FRAGPIPESDRF

Write the FragPipe SDRF template (`sdrf.tsv`) for a run.

## Description

FRAGPIPESDRF ports the FragPipe 24.0 writer triggered by `workflow.misc.save-sdrf=true` (`FragpipeRun.java`, `TabWorkflow.sdrfSave`, `SDRFtable.java`, TabMsfragger SDRF helpers).
The organism is guessed from UniProt `OS=` headers of the searched FASTA, and modification names come from `unimod.obo`.

## Inputs

| Channel                                                                       | Type  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ----------------------------------------------------------------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, workflow_file, manifest, searched_names, fasta, unimod, logs, tmt_label | tuple | `workflow_file`: the saved `fragpipe.workflow` (all `msfragger.*` keys, `workflow.misc.sdrf-type`); `manifest`: `fragpipe-files.fp-manifest`; `searched_names`: TSV mapping each manifest file name to the file MSFragger searched; `fasta`: searched database; `unimod`: `unimod.obo` of the FragPipe tools folder; `logs`: MSFragger log(s) (calibrated fragment tolerance, instrument); `tmt_label`: TMT-Integrator QuantLabel name (e.g. `TMT-16`) or an empty string |

## Outputs

| Channel         | Type                    | Description    |
| --------------- | ----------------------- | -------------- |
| sdrf            | tuple(val, path)        | `sdrf.tsv`     |
| versions_python | tuple (topic: versions) | Python version |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.args                              | ''      | Additional CLI arguments              |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in `workflows/fragpipe.nf`, `workflows/tmt_labelcheck.nf` and `subworkflows/local/fragpipe_search/main.nf`.

## References

- [SDRF-Proteomics](https://github.com/bigbio/proteomics-sample-metadata)

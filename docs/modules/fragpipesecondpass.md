# FRAGPIPESECONDPASS

FragPipe's second-pass workflow and manifest after a Write sub mzML run.

## Description

When `tab-run.write_sub_mzml=true`, FragPipe's finalizer writes `fragpipe-second-pass.workflow` (Write sub mzML off, Percolator on, PeptideProphet off, Print decoys off, no mass calibration, MSFragger's optimized values parsed from the console text) and `fragpipe-files-second-pass.fp-manifest` (naming `<stem>_sub.mzML`), see `FragpipeRun.java:589-700`.
`bin/fragpipe_second_pass.py` writes both from the staged files.

## Inputs

| Channel                                      | Type  | Description                                                                                                                                                                                                                                                 |
| -------------------------------------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, workflow_file, manifest, logs, workdir | tuple | `workflow_file`: the run's `fragpipe.workflow`; `manifest`: the run's `fragpipe-files.fp-manifest`; `logs`: MSFragger console output in FragPipe's order; `workdir`: FragPipe's work dir (the published `fragpipe/` folder) where the `_sub.mzML` files are |

## Outputs

| Channel         | Type                    | Description                              |
| --------------- | ----------------------- | ---------------------------------------- |
| workflow        | tuple(val, path)        | `fragpipe-second-pass.workflow`          |
| manifest        | tuple(val, path)        | `fragpipe-files-second-pass.fp-manifest` |
| versions_python | tuple (topic: versions) | Python version                           |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.args                              | ''      | Additional CLI arguments              |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in `workflows/fragpipe.nf` and `workflows/tmt_labelcheck.nf` when Write sub mzML is on.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

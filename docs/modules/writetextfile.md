# WRITETEXTFILE

Write a text file whose body was built outside the tool modules.

## Description

WRITETEXTFILE writes a file body built by Groovy helpers or the workflow parser, so tool modules receive it as a path and only run their tool.
The body is written byte for byte (a body that does not end with a newline gets none).
The process is included under aliases per use (for example `WRITE_FRAGGER_PARAMS`, `WRITE_LABELQUANT_ANNOTATION`).

## Inputs

| Channel          | Type                 | Description                                              |
| ---------------- | -------------------- | -------------------------------------------------------- |
| meta, name, text | tuple(val, val, val) | Metadata map carried to the output, file name, file body |

## Outputs

| Channel            | Type                    | Description       |
| ------------------ | ----------------------- | ----------------- |
| file               | tuple(val, path)        | The written file  |
| versions_coreutils | tuple (topic: versions) | coreutils version |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`quay.io/nf-core/ubuntu:20.04`

## Usage

Used in the `fragpipe_search`, `fragpipe_quant`, `fragpipe_philosopher_quant` and `fragpipe_speclib` subworkflows.

## References

- [Nextflow processes](https://www.nextflow.io/docs/latest/process.html)

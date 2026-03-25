# FRAGPIPE_HEADLESS

Run FragPipe in headless mode as a single all-in-one process.

## Description

FRAGPIPE_HEADLESS runs FragPipe in headless mode, producing results identical to the GUI.
Unlike the FRAGPIPE module, this one takes raw manifest content as a val input and stages raw files into a subdirectory, generates the manifest from provided content, and updates the workflow file for the container environment by stripping developer-machine tool paths.
Used in the main `msproteomics.nf` workflow as an alternative to the modular FRAGPIPE_WF subworkflow.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| raw_files | path | All raw data files (`.d`, `.raw`, `.mzML`) staged into `raw_files/` subdirectory |
| database | path | FASTA database file for protein identification |
| workflow_file | path | FragPipe `.workflow` configuration file |
| manifest_content | val(string) | Tab-separated manifest content (filename, experiment, bioreplicate, data_type) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| all_results | path | All FragPipe output files (`results/**`) |
| combined_protein | path | Combined protein report across all experiments (`combined_protein.tsv`, optional) |
| combined_peptide | path | Combined peptide report across all experiments (`combined_peptide.tsv`, optional) |
| combined_ion | path | Combined ion report across all experiments (`combined_ion.tsv`, optional) |
| versions | path | Software versions (`versions.yml`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| (memory) | 90% of task.memory | Java heap size computed automatically from allocated task memory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the main workflow (`workflows/msproteomics.nf`) as a standalone all-in-one FragPipe execution.
Provides an alternative to the modular subworkflow approach when the full FragPipe pipeline should run as a single process.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [FragPipe tutorial](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html)

# SAINTEXPRESS

Scores protein-protein interactions from AP-MS data using SAINTexpress.

## Description

SAINTexpress is a statistical tool for evaluating the significance of prey proteins identified in affinity purification-mass spectrometry (AP-MS) experiments.
It supports both spectral count (`spc`) and intensity (`int`) scoring modes.
This module discovers the SAINTexpress binary from the FragPipe tools directory and runs it on interaction, bait, and prey definition files.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| inter_file | path | Interaction file containing prey-bait spectral counts or intensities (`inter.dat`) |
| bait_file | path | Bait definition file mapping samples to baits and controls (`bait.dat`) |
| prey_file | path | Prey definition file with protein lengths (`prey.dat`) |
| mode | val(string) | Scoring mode: `"spc"` for spectral count or `"int"` for intensity |
| config_cli | val(string) | SAINTexpress parameters as CLI string (e.g., `-R 2 -L 4`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing SAINTexpress results and log files |
| results_list | tuple(val, path) | `${prefix}/list.txt` with interaction scores |
| versions_saintexpress | topic: versions | SAINTexpress version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended after `config_cli` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory for binary discovery |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`) for AP-MS protein-protein interaction scoring.

## References

- [SAINTexpress homepage](https://saint-apms.sourceforge.net/Main.html)
- [Publication](https://doi.org/10.1016/j.jprot.2013.10.023)

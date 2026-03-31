# SKYLINE

Creates Skyline documents from FragPipe results for targeted quantification.

## Description

Skyline is an open-source application for building targeted mass spectrometry methods and analyzing the resulting data.
This module creates Skyline documents from FragPipe search results via the FragPipe Java wrapper (`org.nesvilab.fragpipe.tools.skyline.Skyline`), importing spectral libraries and mass spectrometry data.
Memory is auto-calculated from `task.memory`.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| results_dir | path | FragPipe output directory containing psm.tsv and library files |
| mzml_files | path (stageAs `spectra/*`) | Mass spectrometry files for import (`*.mzML`, `*.mzXML`) |
| speclib | path | Spectral library file (`*.speclib`, `*.sptxt`, `*.blib`) |
| config_cli | val(string) | Skyline parameters as CLI string (modsMode, tolerances, etc.) |
| skyline_path | val(string) | Path to Skyline executable (shared resource) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing Skyline output files |
| skyline_document | tuple(val, path) | `${prefix}/skyline_files/fragpipe.sky` Skyline document (optional) |
| reports | tuple(val, path) | `${prefix}/*.csv` quantification report files (optional) |
| versions_skyline | topic: versions | Skyline version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for the Skyline Java wrapper |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory for JAR discovery |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`) for creating Skyline documents from FragPipe search results.

## References

- [Skyline homepage](https://skyline.ms/)
- [Skyline tutorials](https://skyline.ms/wiki/home/software/Skyline/page.view?name=tutorials)
- [Publication](https://doi.org/10.1093/bioinformatics/btq054)

# MBG

Mass-based glycoproteomics matching for glycan identification.

## Description

MBG (Mass-Based Glycoproteomics) performs glycan identification by analyzing PSM data against glycan databases, matching observed mass shifts to known glycan compositions.
It supports both standard and TMT-based glycoproteomics workflows, discovers the MBG JAR and dependencies from the FragPipe tools directory, and supports Bruker and Thermo native libraries.
Used within the FRAGPIPE_GLYCO subworkflow for glycoproteomics analysis.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| psm_file | path | PSM file from PHILOSOPHER_FILTER (`psm.tsv`) |
| manifest_file | path | FragPipe manifest file listing input files (`*.fp-manifest`) |
| config_cli | val(string) | MBG parameters as CLI string (e.g., `--maxq 0.01 --mztol 10`) |
| residue_db | path | Glycan residues database file (shared resource) |
| glycan_mod_db | path | Glycan modifications database file (shared resource) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing MBG results and log |
| results | tuple(meta, path) | MBG matching results (`*_mbg*.tsv`, optional) |
| glycan_results | tuple(meta, path) | Glycan-level results (`*_glycan*.tsv`, optional) |
| versions_mbg | tuple | MBG software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the MBG Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_GLYCO` subworkflow for glycoproteomics analysis within the modular FragPipe pipeline.
Results are emitted as `mbg_results` from the `FRAGPIPE_WF` subworkflow (`subworkflows/local/fragpipe/main.nf`).

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [Glycoproteomics tutorial](https://fragpipe.nesvilab.org/docs/tutorial_glyco.html)
- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

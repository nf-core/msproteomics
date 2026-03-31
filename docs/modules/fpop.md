# FPOP

Perform Fast Photochemical Oxidation of Proteins (FPOP) analysis.

## Description

FPOP analysis calculates oxidation levels per residue/region, compares FPOP vs control samples, and supports both label-free (LFQ) and TMT-based workflows.
The analysis is performed using the `FragPipe_FPOP_Analysis.py` Python script bundled with FragPipe, with configuration provided as bash-sourceable `KEY='value'` lines.
Used within the FRAGPIPE_EXPORT subworkflow for structural proteomics experiments.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| input_file | path | Peptide quantification file (`combined_modified_peptide.tsv` for LFQ, or TMT abundance file) |
| secondary_file | path | Optional secondary TMT file (use `NO_FILE` placeholder if not applicable) |
| config_cli | val(string) | Bash-sourceable config: REGION_SIZE, CONTROL_LABEL, FPOP_LABEL, SUBTRACT_CONTROL, IS_TMT |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing FPOP analysis results and log |
| results | tuple(meta, path) | FPOP analysis results in TSV format (`*_fpop*.tsv`, optional) |
| results_csv | tuple(meta, path) | FPOP analysis results in CSV format (`*_fpop*.csv`, optional) |
| versions_fpop | tuple | FPOP software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the FPOP Python command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`) for FPOP structural proteomics analysis.
Enabled conditionally via `shouldRunTool(configs, 'fpop')`.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [FPOP tutorial](https://fragpipe.nesvilab.org/docs/tutorial_fpop.html)
- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

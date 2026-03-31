# IONQUANT

Label-free and isobaric quantification with match-between-runs and MaxLFQ.

## Description

IonQuant is a fast and comprehensive tool for MS1 quantification of label-free and isobaric labeled proteomics data, providing accurate feature detection, match-between-runs (MBR), and MaxLFQ normalization.
It operates as an aggregate step across all samples, generating a filelist pointing to per-sample psm.tsv files and spectra directories, and supports TMT isobaric quantification via annotation files.
Used within the FRAGPIPE_QUANT subworkflow and the TMT Label Check workflow for quantification.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| psm_dirs | path | PHILOSOPHER_FILTER output directories (named by sample, containing `psm.tsv`) |
| spec_files | path | Spectral files (`*.mzML` or `*.raw`) staged in `spectra/` subdirectory |
| annotation_file | path | TMT annotation file (2-column or 6-column TSV; pass `[]` if not applicable) |
| config_cli | val(string) | IonQuant parameters as CLI string (e.g., `--mbr 1 --maxlfq 1`) |
| modmasses | val(string) | Comma-separated modification masses for `--modlist` (empty = skip) |
| ionquant_dir | path | Unzipped IonQuant tool directory (optional, pass `[]` when not using) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing all IonQuant output files |
| ions | tuple(meta, path) | Ion-level quantification (`ion.tsv`, optional) |
| combined_ions | tuple(meta, path) | Combined ion-level quantification across samples (`combined_ion.tsv`, optional) |
| combined_protein | tuple(meta, path) | Combined protein-level quantification (`combined_protein.tsv`, optional) |
| combined_peptide | tuple(meta, path) | Combined peptide-level quantification (`combined_peptide.tsv`, optional) |
| combined_modified_peptide | tuple(meta, path) | Combined modified-peptide-level quantification (`combined_modified_peptide.tsv`, optional) |
| combined_site | tuple(meta, path) | Combined site-level quantification (`combined_site_*.tsv`, optional) |
| all_tsv | tuple(meta, path) | All TSV output files |
| license_agreement | path | License agreement marker file |
| versions_ionquant | tuple | IonQuant software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the IonQuant Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` or the process exits with an error |

## Container

Inherited from pipeline configuration (requires FragPipe container).

## Usage

Used in the `FRAGPIPE_QUANT` subworkflow (`subworkflows/local/fragpipe_quant/main.nf`) as `IONQUANT`, `IONQUANT_MS1`, and `IONQUANT_ISOBARIC` for different quantification modes.
Also referenced in the TMT Label Check workflow (`subworkflows/local/tmt_labelcheck/main.nf`).
Includes a guard that skips execution when both `--perform-isoquant 0` and `--perform-ms1quant 0` are set.

## References

- [IonQuant homepage](https://ionquant.nesvilab.org/)
- [IonQuant GitHub](https://github.com/Nesvilab/IonQuant)
- [IonQuant Wiki](https://github.com/Nesvilab/IonQuant/wiki)
- [Publication: doi:10.1074/mcp.TIR120.002048](https://doi.org/10.1074/mcp.TIR120.002048)

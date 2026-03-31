# MSFRAGGER

Ultrafast peptide identification using MSFragger database search engine.

## Description

MSFragger is an ultrafast database search tool for peptide identification in mass spectrometry-based proteomics.
It uses a fragment-ion indexing method to achieve speeds 10-100x faster than conventional search engines.
This module supports both params-file mode (FragPipe calling convention) and CLI mode, handles split-database search via `--partial` flag, and auto-discovers MSFragger JAR and native libraries (Bruker, Thermo).

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (may include `chunk_id` for split mode) |
| mzml_files | path | Mass spectrometry data files (`*.mzML`) |
| fasta | path | Protein sequence database (`*.fasta`) |
| params_file | path | MSFragger params file (native fragger.params format). Pass `[]` for CLI mode |
| pepindex | path | Prebuilt pepindex files from MSFRAGGER_INDEX. Pass `[]` when not using |
| msfragger_dir | path | Unzipped MSFragger tool directory (optional). Pass `[]` when not using |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| pepxml | tuple(val, path) | PepXML search results (`*.pepXML`, optional) |
| pin | tuple(val, path) | Percolator input files (`*.pin`, optional) |
| tsv | tuple(val, path) | TSV outputs such as mass calibration results (optional) |
| calibrated_mzml | tuple(val, path) | Calibrated mzML files (when write_calibrated_mzml=true, optional) |
| uncalibrated_mzml | tuple(val, path) | Uncalibrated mzML for .d/.raw inputs (optional) |
| results_dir | tuple(val, path) | Results directory containing all output files |
| log | tuple(val, path) | Log files (optional) |
| license_agreement | path | License agreement file (`I_AGREE_FRAGPIPE_LICENSE_AGREEMENT`) |
| versions_msfragger | tuple (topic: versions) | MSFragger version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the MSFragger Java command |
| ext.calibrate_mass | (not set) | Override calibrate_mass in params file (only when explicitly set) |
| ext.write_calibrated_mzml | (not set) | Override write_calibrated_mzml (only when explicitly set) |
| ext.java_xmx | (auto) | Override Java heap size |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` to run |

## Container

Not defined in module; inherited from workflow/config (typically `docker.io/fcyucn/fragpipe:24.0`).

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`).
In split-database mode (`meta.chunk_id` set), outputs are prefixed with `chunk_<id>_` to avoid collisions.
Requires explicit license agreement via `ext.agree_fragpipe_license_agreement = true`.

## References

- [MSFragger homepage](https://msfragger.nesvilab.org/)
- [MSFragger Wiki](https://github.com/Nesvilab/MSFragger/wiki)
- [Publication: Kong et al., Nature Methods, 2017](https://doi.org/10.1038/nmeth.4256)

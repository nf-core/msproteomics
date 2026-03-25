# MSFRAGGER_CALIBRATE

Mass calibration and parameter optimization for split-database search using MSFragger.

## Description

MSFRAGGER_CALIBRATE replicates the `calibrate()` function from FragPipe's `msfragger_pep_split.py`.
It runs MSFragger with the `--split1` flag on ALL samples to perform mass calibration, then extracts optimized parameters (fragment_mass_tolerance, etc.) via `parse_calibration_output.py`.
This is an AGGREGATE process that runs once with all samples collected together, producing calibrated spectra and an updated fragger.params file for downstream split-database search.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| fasta | path | Protein sequence database (`*.fasta`) |
| mzml_files | path | ALL mass spectrometry data files collected (`*.mzML`) |
| params_file | path | MSFragger params file (native fragger.params format) |
| msfragger_dir | path | Unzipped MSFragger tool directory (optional). Pass `[]` when not using |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| calibrated_spectra | path | `.mzBIN_calibrated` files (or original copies if not calibrated) |
| params | path | Updated `calibrated_fragger.params` with optimized tolerances, `check_spectral_files=0`, `calibrate_mass=0` |
| license_agreement | path | License agreement file (`I_AGREE_FRAGPIPE_LICENSE_AGREEMENT`) |
| versions_msfragger | tuple (topic: versions) | MSFragger version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the MSFragger `--split1` command |
| ext.java_xmx | (auto) | Override Java heap size |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` to run |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

Not defined in module; inherited from workflow/config (typically `docker.io/fcyucn/fragpipe:24.0`).

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) when split-database search is enabled.
Steps: (1) sort FASTA proteins lexicographically, (2) create calibration params, (3) run MSFragger `--split1`, (4) parse output via `parse_calibration_output.py`.

## References

- [MSFragger homepage](https://msfragger.nesvilab.org/)
- [MSFragger Wiki](https://github.com/Nesvilab/MSFragger/wiki)
- [Publication: Kong et al., Nature Methods, 2017](https://doi.org/10.1038/nmeth.4256)

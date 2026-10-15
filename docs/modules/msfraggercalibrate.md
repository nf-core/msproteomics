# MSFRAGGERCALIBRATE

Mass calibration and parameter optimization for split-database search using MSFragger.

## Description

MSFRAGGERCALIBRATE replicates `calibrate()` of FragPipe 24.0 `tools/msfragger_pep_split.py` (lines 498-611) for one spectra file and its data type's params set.
FragPipe runs the split script once per params set; the per-file patched FragPipe used for golden comparison runs it per file, so this process runs once per file whose params have `calibrate_mass` 1 or 2 (files with `calibrate_mass = 0` are searched without calibration, `msfragger_pep_split.py:630-632`).
It runs MSFragger `--split1` against the sorted full FASTA, then `parse_calibration_output.py` applies the optimized values (only for `calibrate_mass = 2`) and `check_spectral_files = 0` to the file's params, keeping every other line (including `calibrate_mass`) as FragPipe does.

## Inputs

| Channel                    | Type                   | Description                                                                             |
| -------------------------- | ---------------------- | --------------------------------------------------------------------------------------- |
| meta, spectra, params_file | tuple(val, path, path) | One spectra file and its MSFragger params (`fragger.params`, `fragger_dia.params`, ...) |
| fasta                      | path                   | Protein sequence database (`*.fasta`)                                                   |
| msfragger_dir              | path                   | Unzipped MSFragger tool directory (optional). Pass `[]` when not using                  |

## Outputs

| Channel            | Type                    | Description                                                                                         |
| ------------------ | ----------------------- | --------------------------------------------------------------------------------------------------- |
| calibrated         | tuple(val, path, path)  | `.mzBIN_calibrated` (or a copy of the input) and the calibrated params (`calibrated/<params name>`) |
| uncalibrated_mzml  | tuple(val, path)        | `*_uncalibrated.mzML` MSFragger writes for `.raw` / `.d` input (optional)                           |
| log                | tuple(val, path)        | MSFragger console output (`msfragger_calibrate.log`), read by FRAGPIPESDRF                          |
| license_agreement  | path                    | License agreement file (`I_AGREE_FRAGPIPE_LICENSE_AGREEMENT`)                                       |
| versions_msfragger | tuple (topic: versions) | MSFragger version                                                                                   |
| versions_fragpipe  | tuple (topic: versions) | FragPipe version                                                                                    |

## Parameters

| Parameter                            | Default        | Description                                                           |
| ------------------------------------ | -------------- | --------------------------------------------------------------------- |
| ext.args                             | ''             | Additional CLI arguments appended to the MSFragger `--split1` command |
| ext.java_xmx                         | (auto)         | Override Java heap size                                               |
| ext.agree_fragpipe_license_agreement | (required)     | Must be `true` to run                                                 |
| ext.fragpipe_tools_dir               | (default path) | Override FragPipe tools directory                                     |

## Container

Not defined in module; inherited from workflow/config (the licensed FragPipe image given with `--fragpipe_container`).

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) when split-database search is enabled.
Steps: (1) sort FASTA proteins lexicographically, (2) create calibration params, (3) run MSFragger `--split1`, (4) parse output via `parse_calibration_output.py`.

## References

- [MSFragger homepage](https://msfragger.nesvilab.org/)
- [MSFragger Wiki](https://github.com/Nesvilab/MSFragger/wiki)
- [Publication: Kong et al., Nature Methods, 2017](https://doi.org/10.1038/nmeth.4256)

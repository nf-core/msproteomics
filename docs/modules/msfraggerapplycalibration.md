# MSFRAGGERAPPLYCALIBRATION

Apply MSFragger's calibration output to the search params and spectra of a split-database search.

## Description

MSFRAGGERAPPLYCALIBRATION runs `parse_calibration_output.py`, the port of `calibrate()` in FragPipe 24.0 `tools/msfragger_pep_split.py` (lines 524-611).
The optimized values MSFragger printed are applied to the search params (`calibrate_mass = 2` only), `check_spectral_files` is set to 0, and the spectra file is replaced by its `.mzBIN_calibrated` (or copied when MSFragger wrote none).

## Inputs

| Channel                                                     | Type  | Description                                                                                                                                                        |
| ----------------------------------------------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| meta, spectra, calibrate_log, calibrated_mzbin, params_file | tuple | The spectra file MSFRAGGERCALIBRATE calibrated, its console output, the `<stem>.mzBIN_calibrated` it wrote (`[]` when none) and the file's MSFragger search params |

## Outputs

| Channel         | Type                    | Description                                                                                               |
| --------------- | ----------------------- | --------------------------------------------------------------------------------------------------------- |
| calibrated      | tuple(val, path, path)  | Calibrated spectra (`calibrated_spectra/*`) and the params after calibration (`calibrated/<params name>`) |
| versions_python | tuple (topic: versions) | Python version                                                                                            |

## Parameters

| Parameter                             | Default | Description                           |
| ------------------------------------- | ------- | ------------------------------------- |
| ext.args                              | ''      | Additional CLI arguments              |
| ext.singularity_pull_docker_container | false   | Use the Docker image with Singularity |
| ext.when                              | (none)  | Skip the process when false           |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) after MSFRAGGERCALIBRATE.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

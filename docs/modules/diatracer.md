# DIATRACER

Convert Bruker `.d` files to mzML for DIA analysis.

## Description

diaTracer processes Bruker timsTOF `.d` files for DIA (Data Independent Acquisition) analysis, converting them to pseudo-DDA mzML format.
It is part of the MSFragger/FragPipe suite and requires a FragPipe license agreement (`ext.agree_fragpipe_license_agreement = true`).
Used within the FRAGPIPE_CONVERT subworkflow for Bruker DIA file conversion.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| d_file | path | Bruker `.d` directory containing raw mass spectrometry data |
| config_cli | val(string) | diaTracer parameters as CLI string |
| diatracer_dir | path | Unzipped diaTracer tool directory (optional, pass `[]` when not using) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing diaTracer results and log |
| mzml | tuple(meta, path) | Converted mzML file (`*_diatracer.mzML`) |
| license_agreement | path | License agreement marker file |
| versions_diatracer | tuple | diaTracer software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the diaTracer Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` or the process exits with an error |

## Container

Inherited from pipeline configuration (requires FragPipe container).

## Usage

Used in the `FRAGPIPE_CONVERT` subworkflow (`subworkflows/local/fragpipe_convert/main.nf`) to convert Bruker `.d` files for DIA analysis.
Enabled conditionally via `shouldRunTool(configs, 'diatracer')`.

## References

- [MSFragger GitHub](https://github.com/Nesvilab/MSFragger)
- [MSFragger Wiki](https://github.com/Nesvilab/MSFragger/wiki)

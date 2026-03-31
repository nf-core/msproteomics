# OPENMS_PEAK_PICKER

Centroid peak picking using OpenMS PeakPickerHiRes.

## Description

OPENMS_PEAK_PICKER converts profile-mode mass spectra to centroid mode using the PeakPickerHiRes algorithm from OpenMS.
Centroided spectra are required by most downstream proteomics analysis tools (database search, quantification, etc.).
This module supports in-memory or low-memory processing modes and configurable MS levels for peak picking.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (uses `meta.mzml_id`) |
| mzml_file | path | Input mzML file in profile mode (`*.mzML`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| mzmls_picked | tuple(val, path) | Centroided mzML file (`*.mzML`) |
| versions | path | PeakPickerHiRes version (`versions.yml`) |
| log | path | Log file (`*.log`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to `PeakPickerHiRes` command |
| params.peakpicking_inmemory | - | Processing mode: `inmemory` or `lowmemory` |
| params.peakpicking_ms_levels | - | MS levels for peak picking (sets `-algorithm:ms_levels`) |
| params.pp_debug | - | Debug level (sets `-debug`) |

## Container

`ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14`

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) to convert profile-mode spectra to centroided format before database search and quantification.

## References

- [PeakPickerHiRes Documentation](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_PeakPickerHiRes.html)
- [OpenMS homepage](https://www.openms.de/)

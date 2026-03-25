# MZML_STATISTICS

Computes statistics from mzML files using quantms-utils.

## Description

MZML_STATISTICS produces MS-level information parquet files and optionally generates MS2-level statistics and feature detection results from mass spectrometry data.
It uses the `quantmsutilsc mzmlstats` command from the quantms-utils package.
Behavior is controlled by `params.id_only` (skips MS2 stats when true) and `params.mzml_features` (enables feature detection).

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (uses `meta.mzml_id`) |
| ms_file | path | Mass spectrometry file in mzML format (`*.mzML`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| ms_statistics | path | MS-level statistics (`*_ms_info.parquet`) |
| ms2_statistics | tuple(val, path) | MS2-level statistics (`*_ms2_info.parquet`, optional) |
| feature_statistics | path | Feature detection results (`*_feature_info.parquet`, optional) |
| versions | path | quantms-utils version (`versions.yml`) |
| log | path | Log file (`*.log`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to `quantmsutilsc mzmlstats` command |
| params.id_only | - | When true, skips MS2-level statistics (`--ms2_file` not added) |
| params.mzml_features | - | When true, enables feature detection (`--feature_detection` flag) |

## Container

`biocontainers/quantms-utils:0.0.24--pyh7e72e81_0`

## Usage

Available as a local module for computing mzML-level QC statistics.
Not currently included in any active workflow but can be used for quality control reporting.

## References

- [quantms GitHub](https://github.com/bigbio/quantms)
- [quantms Documentation](https://github.com/bigbio/quantms/tree/readthedocs)

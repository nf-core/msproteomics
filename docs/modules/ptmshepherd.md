# PTMSHEPHERD

Post-translational modification analysis and mass shift characterization using PTM-Shepherd.

## Description

PTM-Shepherd is a bioinformatics tool for characterization of post-translational and chemical modifications from open search results.
It performs mass shift profiling, modification summarization, diagnostic ion mining, PTM localization, and glycoprofile analysis from PSM-level results.
PTM-Shepherd is part of the FragPipe suite and uses multiple JAR dependencies including batmass-io, commons-math3, hipparchus-core, hipparchus-stat, and optionally IonQuant.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| psm_dirs | path | Directories containing `psm.tsv` files from PHILOSOPHER_FILTER |
| protxml | path | Combined protein inference file from PROTEINPROPHET (`*.prot.xml`) |
| mzml_files | path | Mass spectrometry files (`*.mzML`, staged in `spectra/` directory) |
| config_file | path | PTM-Shepherd configuration file (`shepherd.config`) |
| fasta | path | Protein database FASTA (shared resource) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Results directory containing all PTM-Shepherd output |
| global_profile | tuple(val, path) | Global mass shift profile across all samples (`global.profile.tsv`, optional) |
| global_modsummary | tuple(val, path) | Global modification summary with annotations (`global.modsummary.tsv`, optional) |
| diagmine | tuple(val, path) | Diagnostic ion mining results (`*diagmine.tsv`, optional) |
| localization | tuple(val, path) | PTM site localization results (`*localization.tsv`, optional) |
| glycoprofile | tuple(val, path) | Glycoproteomics profile results (`*glycoprofile.tsv`, optional) |
| versions_ptmshepherd | tuple (topic: versions) | PTM-Shepherd version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to PTM-Shepherd Java command after config file path |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_PTM` subworkflow (`subworkflows/local/fragpipe_ptm/main.nf`) for post-search PTM characterization.
The module injects `threads`, `database`, and `dataset` lines into `shepherd.config`, overriding any existing values.
Dataset lines are built from per-sample PSM directories and spectra files.
Memory is auto-calculated from `task.memory` (minus 2 GB headroom).

## References

- [PTM-Shepherd GitHub](https://github.com/Nesvilab/PTM-Shepherd)
- [PTM-Shepherd Wiki](https://github.com/Nesvilab/PTM-Shepherd/wiki)
- [Publication: Geiszler et al., J. Proteome Research, 2021](https://doi.org/10.1021/acs.jproteome.0c00688)

# MSFRAGGER_INDEX

Prebuild MSFragger pepindex files from FASTA and params (digest-only mode).

## Description

MSFRAGGER_INDEX runs MSFragger with a params file and FASTA but no spectra files, which triggers digest-only mode and creates pepindex files alongside the FASTA.
This replicates `CmdMsfraggerDigest.java` from FragPipe: the pepindex files cache the in-silico digest and can be reused by downstream MSFRAGGER search tasks to skip redundant digest computation.
The module sets `calibrate_mass = 0` in the params to prevent calibration during indexing.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata map |
| fasta | path | Protein sequence database (`*.fasta`) |
| params_file | path | MSFragger params file (native fragger.params format) |
| msfragger_dir | path | Unzipped MSFragger tool directory (optional). Pass `[]` when not using |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| indexed_fasta | tuple(val, path, path) | FASTA + co-located pepindex files (must stay together for MSFragger reuse) |
| license_agreement | path | License agreement file (`I_AGREE_FRAGPIPE_LICENSE_AGREEMENT`) |
| versions_msfragger | tuple (topic: versions) | MSFragger version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the MSFragger digest command |
| ext.java_xmx | (auto) | Override Java heap size |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` to run |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

Not defined in module; inherited from workflow/config (typically `docker.io/fcyucn/fragpipe:24.0`).

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) to prebuild peptide indexes that are reused across multiple search tasks, avoiding redundant digest computation.
Pepindex depends on digest parameters only: enzyme, missed cleavages, peptide length/mass range, and variable modifications.

## References

- [MSFragger homepage](https://msfragger.nesvilab.org/)
- [MSFragger Wiki](https://github.com/Nesvilab/MSFragger/wiki)
- [Publication: Kong et al., Nature Methods, 2017](https://doi.org/10.1038/nmeth.4256)

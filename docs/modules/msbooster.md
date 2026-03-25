# MSBOOSTER

Deep learning-based PSM rescoring using predicted retention time and spectra.

## Description

MSBooster enhances peptide identification by using deep learning predictions of retention time and MS2 spectra from DIA-NN.
It rescores PSMs by comparing observed vs predicted properties, adding features to Percolator input (PIN) files for improved sensitivity.
MSBooster is part of the FragPipe suite and fits in the pipeline between MSFragger database search and Percolator validation.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| pin_files | path | PIN files from MSFragger search (`*.pin`) |
| mzml_files | path | Mass spectrometry data files (`*.mzML`) |
| params_file | path | Full MSBooster params file (key=value format). Pass `[]` when not using |
| fragger_params | path | MSFragger params file for modification definitions. Pass `[]` when not using |
| has_ion_mobility | val(boolean) | Whether input data has ion mobility (e.g., Bruker timsTOF .d files) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Results directory containing rescored PIN files and logs |
| pin_edited | tuple(val, path) | Rescored PIN files with additional features (renamed from `*_edited.pin`) |
| versions_msbooster | tuple (topic: versions) | MSBooster version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to `java ... MainClass --paramsList` command |
| ext.fragpipe_tools_dir | '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools' | Override FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`).
Supports two modes: params-file mode (matching FragPipe calling convention) and CLI mode (backward compatible).
In params-file mode, runtime parameters (numThreads, DiaNN, unimodObo, fragger, mzmlDirectory, pinPepXMLDirectory) are injected/overridden.

## References

- [MSBooster GitHub](https://github.com/Nesvilab/MSBooster)
- [MSBooster Wiki](https://github.com/Nesvilab/MSBooster/wiki)

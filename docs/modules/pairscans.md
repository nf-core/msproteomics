# PAIRSCANS

Pair parent and child scans of one LC-MS file for O-Pair.

## Description

PAIRSCANS runs FragPipe's `PairScans` utility (`CmdPairScans.java`, `util/PairScans.java`, class `org.nesvilab.fragpipe.util.PairScans`) on one LC-MS file.
The `.pairs` file is written next to the spectra file as PairScans does; PairScans writes nothing for an unknown activation type or format (`PairScans.java:80-88`), and O-Pair then treats the file as single scan type, as in FragPipe.

## Inputs

| Channel                       | Type                  | Description                                                                                                                                                   |
| ----------------------------- | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, spectra, pairscans_args | tuple(val, path, val) | `spectra`: the LC-MS file FragPipe lists (mzML, mzXML or Thermo `.raw`); `pairscans_args`: `[activation1, activation2, reverse_scan_order, single_scan_type]` |

## Outputs

| Channel           | Type                    | Description                           |
| ----------------- | ----------------------- | ------------------------------------- |
| pairs             | tuple(val, path)        | `spectra/<basename>.pairs` (optional) |
| log               | tuple(val, path)        | `<basename>.pairscans.log`            |
| versions_fragpipe | tuple (topic: versions) | FragPipe version                      |

## Parameters

| Parameter  | Default | Description                                 |
| ---------- | ------- | ------------------------------------------- |
| ext.args   | ''      | Additional CLI arguments                    |
| ext.prefix | meta.id | Override the output directory / file prefix |
| ext.when   | (none)  | Skip the process when false                 |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in the `FRAGPIPE_GLYCO` subworkflow (`subworkflows/local/fragpipe_glyco/main.nf`) before OPAIR.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)

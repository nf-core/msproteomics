# CHECKCENTROID

FragPipe's centroid check (`CmdCheckCentroid`, `org.nesvilab.fragpipe.util.CheckCentroid`).

## Description

FragPipe runs the check once, before DIA-Umpire, diaTracer and MSFragger (`FragpipeRun.java:1333-1340`), on the first LC-MS file: the first file of the alphabetically first group, without DIA-Quant files (`checkCentroidFile`).
An mzML or mzXML file with a profile-mode scan fails the run; other formats pass (`CheckCentroid.java:41-58`).
The search waits for the check.
argv: `java -Xmx<ram>G -cp <tools>/../lib/fragpipe-24.0.jar:<tools>/batmass-io-1.36.6.jar org.nesvilab.fragpipe.util.CheckCentroid <file> <threads>`.

## Inputs

| Channel | Type     | Description                     |
| ------- | -------- | ------------------------------- |
| meta    | val(map) | Metadata of the checked file    |
| lcms    | path     | The first LC-MS file of the run |

## Outputs

| Channel           | Type             | Description                                        |
| ----------------- | ---------------- | -------------------------------------------------- |
| log               | tuple(val, path) | `check_centroid.log`: command line and tool output |
| versions_fragpipe | topic: versions  | FragPipe version                                   |

## Parameters

| Parameter              | Default                                           | Description              |
| ---------------------- | ------------------------------------------------- | ------------------------ |
| ext.args               | ''                                                | Extra arguments          |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

The licensed FragPipe image given with `--fragpipe_container`.

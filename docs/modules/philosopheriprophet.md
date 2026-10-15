# PHILOSOPHERIPROPHET

Philosopher iProphet over all validated pepXMLs.

## Description

PHILOSOPHERIPROPHET mirrors FragPipe `CmdIprophet.java`: it runs `philosopher iprophet` over the validated pepXMLs of every group (`sharedPepxmlFiles`, `CmdIprophet.java:53`) and writes `combined.pep.xml`.

## Inputs

| Channel                                        | Type  | Description                                                                                                                                                                                                                                                                            |
| ---------------------------------------------- | ----- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, pepxml_files, pepxml_paths, iprophet_cli | tuple | `pepxml_files`: validated pepXMLs of every group; `pepxml_paths`: each pepXML's FragPipe location relative to the work dir (`<group>/<file>`, or `<file>` for the blank group) in FragPipe order; `iprophet_cli`: `--decoy <tag> --nonsp --output combined` (`CmdIprophet.java:44-50`) |

## Outputs

| Channel              | Type                    | Description                         |
| -------------------- | ----------------------- | ----------------------------------- |
| results_dir          | tuple(val, path)        | iProphet outputs and `iprophet.log` |
| pepxml               | tuple(val, path)        | `combined.pep.xml`                  |
| versions_philosopher | tuple (topic: versions) | Philosopher version                 |
| versions_fragpipe    | tuple (topic: versions) | FragPipe version                    |

## Parameters

| Parameter  | Default | Description                                 |
| ---------- | ------- | ------------------------------------------- |
| ext.args   | ''      | Additional CLI arguments                    |
| ext.prefix | meta.id | Override the output directory / file prefix |
| ext.when   | (none)  | Skip the process when false                 |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in the `FRAGPIPE_INFERENCE` subworkflow (`subworkflows/local/fragpipe_inference/main.nf`); the module lives in `modules/local/philosopheriprophet`.

## References

- [Philosopher](https://github.com/Nesvilab/philosopher)

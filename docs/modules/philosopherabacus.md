# PHILOSOPHERABACUS

Philosopher abacus multi-experiment report.

## Description

PHILOSOPHERABACUS mirrors FragPipe `CmdPhilosopherAbacus.java:73-146`.
It runs in the protXML's directory (FragPipe's work dir) on the group directories and writes `combined_protein.tsv`, `combined_peptide.tsv` and the REPRINT tables.
Abacus reads the database of the first group (Philosopher `lib/aba/pro.go`).

## Inputs

| Channel                                                                         | Type  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| meta, group_workspaces, group_names, protxml, combined_pepxml, abacus_cli, plex | tuple | `group_workspaces`: each group's `.meta` workspace after filter/report, in group name order; `group_names`: the group directory names in the same order; `protxml`: `combined.prot.xml` (`[]` when ProteinProphet did not run); `combined_pepxml`: iProphet `combined.pep.xml` for `--peptide`, else `[]`; `abacus_cli`: `<--picked/--uniqueonly> --razor --reprint --tag <tag> [--protein] [--peptide]`; `plex`: TMT plex for `--labels --plex <n>` (0 otherwise, `CmdPhilosopherAbacus.java:125-129`) |

## Outputs

| Channel              | Type                    | Description                       |
| -------------------- | ----------------------- | --------------------------------- |
| results_dir          | tuple(val, path)        | Result directory                  |
| reprint_spc          | tuple(val, path)        | `reprint.spc.tsv` (optional)      |
| reprint_int          | tuple(val, path)        | `reprint.int.tsv` (optional)      |
| combined_protein     | tuple(val, path)        | `combined_protein.tsv` (optional) |
| combined_peptide     | tuple(val, path)        | `combined_peptide.tsv` (optional) |
| versions_philosopher | tuple (topic: versions) | Philosopher version               |
| versions_fragpipe    | tuple (topic: versions) | FragPipe version                  |

## Parameters

| Parameter  | Default | Description                                 |
| ---------- | ------- | ------------------------------------------- |
| ext.args   | ''      | Additional CLI arguments                    |
| ext.prefix | meta.id | Override the output directory / file prefix |
| ext.when   | (none)  | Skip the process when false                 |

## Container

The licensed FragPipe image given with `--fragpipe_container` (no container is defined in the module).

## Usage

Used in the `FRAGPIPE_INFERENCE` subworkflow (`subworkflows/local/fragpipe_inference/main.nf`).

## References

- [Philosopher](https://github.com/Nesvilab/philosopher)

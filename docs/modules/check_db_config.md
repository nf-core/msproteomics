# CHECKDBCONFIG

Replicates FragPipe's database decoy check (`FragpipeRun.checkDbConfig`).

## Description

FragPipe checks the database before PeptideProphet, Percolator or the Philosopher filter runs.
It counts the FASTA entries whose first header field starts with the workflow decoy tag (`database.decoy-tag`).
In headless mode FragPipe refuses to run when the decoy share is 0, 100 %, or outside 40-60 %.
Databases of 1 GiB or more are not inspected (`TabDatabase.databaseSizeLimit`).
This module applies the same rule and fails with FragPipe's messages.
`FRAGPIPE_DATABASE` runs it on every database: `--database` as given, the `--add_philosopher_decoy` / `--add_philosopher_contaminants` output, and the headless path.

| Condition                  | Error message                          |
| -------------------------- | -------------------------------------- |
| No FASTA entries           | `Database looks to be empty`           |
| 0 % entries with the tag   | `No decoys found in the FASTA file.`   |
| 100 % entries with the tag | `All FASTA entries seem to be decoys.` |
| Below 40 % or above 60 %   | `FASTA file contains <pct> decoys.`    |

## Inputs

| Channel   | Type        | Description                                                                       |
| --------- | ----------- | --------------------------------------------------------------------------------- |
| meta      | val(map)    | Database metadata                                                                 |
| fasta     | path        | Protein database in FASTA format                                                  |
| decoy_tag | val(string) | Decoy protein prefix from the workflow `database.decoy-tag` (via `getDecoyTag()`) |

## Outputs

| Channel         | Type             | Description                                             |
| --------------- | ---------------- | ------------------------------------------------------- |
| fasta           | tuple(val, path) | The input database, unchanged, when it passes the check |
| versions_python | topic: versions  | Python version                                          |

## Parameters

| Parameter | Default | Description                                   |
| --------- | ------- | --------------------------------------------- |
| ext.args  | ''      | Additional arguments for `check_db_config.py` |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## References

- `FragpipeRun.java:2703-2765` (`checkDbConfig`)
- `TabDatabase.java:79, 278-285` (size limit, decoy counting)
- `FastaUtils.java:44-85` (`getDecoysCnt`, `readFasta`)

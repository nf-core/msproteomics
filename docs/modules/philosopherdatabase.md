# PHILOSOPHERDATABASE

FragPipe's database update (`CmdDatabaseUpdate`), run only with `--add_philosopher_decoy` / `--add_philosopher_contaminants`.

## Description

FragPipe headless searches the workflow's `database.db-path` as given, so by default the pipeline does not run this module.
With `--add_philosopher_decoy` / `--add_philosopher_contaminants` it runs FragPipe's "Add decoys" argv, `philosopher database --custom <fasta> [--nodecoys] [--contam --contamprefix]` (`CmdDatabaseUpdate.java:44-53`), inside a Philosopher workspace.
There is no `--prefix`, so Philosopher writes its default `rev_` decoys; FragPipe's decoy check (FragpipeRun.checkDbConfig, in FRAGPIPEPLAN's dry run) then rejects the database when the workflow decoy tag is not `rev_`, as FragPipe does.
With `--custom`, Philosopher adds the cRAP contaminants without a `contam_` prefix (that tag only applies to proteome-ID downloads).

## Inputs

| Channel | Type     | Description               |
| ------- | -------- | ------------------------- |
| meta    | val(map) | Database metadata         |
| fasta   | path     | Target protein FASTA file |

## Outputs

| Channel              | Type             | Description                                                           |
| -------------------- | ---------------- | --------------------------------------------------------------------- |
| fasta                | tuple(val, path) | `${prefix}_philosopher.fasta` database with decoys/contaminants added |
| log                  | tuple(val, path) | `philosopher_database.log`: command line, then the tool output        |
| versions_philosopher | topic: versions  | Philosopher version                                                   |
| versions_fragpipe    | topic: versions  | FragPipe version                                                      |

## Parameters

| Parameter              | Default (`conf/modules.config`)                               | Description                                               |
| ---------------------- | ------------------------------------------------------------- | --------------------------------------------------------- |
| ext.args               | from the two params (`--nodecoys`, `--contam --contamprefix`) | Arguments after `database --custom <fasta>`               |
| ext.prefix             | `${meta.id}`                                                  | Output file prefix                                        |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools`             | FragPipe tools directory for Philosopher binary discovery |

## Container

The licensed FragPipe image given with `--fragpipe_container`.

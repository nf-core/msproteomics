/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_DATABASE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validates the protein database the way FragPipe headless does: the FASTA is searched as
    given and checked for decoys (FragpipeRun.checkDbConfig). Only a FASTA the user asked to
    update (--add_philosopher_decoy / --add_philosopher_contaminants) first gets FragPipe's database update (CmdDatabaseUpdate).

    Modules:
    - PHILOSOPHERDATABASE (ch_fasta only: philosopher database --custom)
    - CHECKDBCONFIG      (every database, FragpipeRun.checkDbConfig)

    Execution: Aggregate (single execution)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHILOSOPHERDATABASE } from '../../../modules/local/philosopherdatabase/main'
include { CHECKDBCONFIG } from '../../../modules/local/checkdbconfig/main'
include { getDecoyTag; shouldRunTool } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_DATABASE {
    take:
    ch_fasta         // channel: [ val(meta), path(fasta) ] - FASTA to update with decoys and contaminants (empty unless --add_philosopher_decoy / --add_philosopher_contaminants)
    ch_prebuilt_db   // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty when ch_fasta is used)
    ch_tool_configs  // channel: val(tool_configs_json) - parsed JSON config (source of database.decoy-tag)

    main:

    // Decoy tag (workflow database.decoy-tag) - single source for every decoy-aware tool
    ch_decoy_tag = ch_tool_configs.map { configs -> getDecoyTag(configs) }

    //
    // FragPipe database update (CmdDatabaseUpdate), only for a FASTA the user asked to update
    //
    PHILOSOPHERDATABASE(ch_fasta)

    // Exactly one of the two sources holds the database.
    ch_database = ch_prebuilt_db
        .concat(PHILOSOPHERDATABASE.out.fasta)
        .first()

    //
    // FragPipe validates the database decoys whenever PeptideProphet, Percolator or the
    // Philosopher filter runs (FragpipeRun.java:1597-1600, 1631-1634, 1704-1707).
    //
    ch_database_branched = ch_database
        .combine(ch_tool_configs)
        .branch { meta, fasta, configs ->
            check: ['peptideprophet', 'percolator', 'filter'].any { tool -> shouldRunTool(configs, tool) }
                return [meta, fasta]
            skip: true
                return [meta, fasta]
        }

    CHECKDBCONFIG(ch_database_branched.check, ch_decoy_tag)

    ch_checked_database = CHECKDBCONFIG.out.fasta
        .mix(ch_database_branched.skip)

    // Value channel for the fasta path that can be broadcast to multiple processes
    ch_fasta_path = ch_checked_database
        .map { _meta, fasta -> fasta }
        .first()

    emit:
    fasta      = ch_checked_database  // channel: [ val(meta), path(fasta) ] - validated database
    fasta_path = ch_fasta_path        // channel: path(fasta) - value channel for broadcasting
    versions   = channel.empty()  // channel: empty - every module reports its version on the 'versions' topic
}

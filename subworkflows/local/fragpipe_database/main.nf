/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_DATABASE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepares protein database with decoys and contaminants using Philosopher.

    Modules:
    - PHILOSOPHER_DATABASE

    Execution: Aggregate (single execution)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHILOSOPHER_DATABASE } from '../../../modules/local/philosopher/database/main'

workflow FRAGPIPE_DATABASE {
    take:
    ch_fasta         // channel: [ val(meta), path(fasta) ] - protein database

    main:

    //
    // Run Philosopher database to add decoys and contaminants
    //
    PHILOSOPHER_DATABASE(ch_fasta)

    // Create a value channel for the fasta path that can be broadcast to multiple processes
    ch_fasta_path = PHILOSOPHER_DATABASE.out.fasta
        .map { _meta, fasta -> fasta }
        .first()

    emit:
    fasta      = PHILOSOPHER_DATABASE.out.fasta  // channel: [ val(meta), path(prepared_fasta) ]
    fasta_path = ch_fasta_path                   // channel: path(fasta) - value channel for broadcasting
}

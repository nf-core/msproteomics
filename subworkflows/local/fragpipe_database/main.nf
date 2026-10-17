/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_DATABASE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The protein database FragPipe searches: the FASTA as given, or, when the user asked to update
    it (--add_philosopher_decoy / --add_philosopher_contaminants), the FASTA after FragPipe's database
    update (CmdDatabaseUpdate). FragPipe runs that update from its Database tab before a run
    (DownloadDbHelper), so it comes before FragPipe's plan, which then checks the decoys itself
    (FragpipeRun.checkDbConfig, FragpipeRun.java:1597-1600, 2703-2760).

    Modules:
    - PHILOSOPHERDATABASE (ch_fasta only: philosopher database --custom)

    Execution: Aggregate (single execution)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHILOSOPHERDATABASE } from '../../../modules/local/philosopherdatabase/main'

workflow FRAGPIPE_DATABASE {
    take:
    ch_fasta         // channel: [ val(meta), path(fasta) ] - FASTA to update with decoys and contaminants (empty unless --add_philosopher_decoy / --add_philosopher_contaminants)
    ch_prebuilt_db   // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty when ch_fasta is used)

    main:
    PHILOSOPHERDATABASE(ch_fasta)

    // Exactly one of the two sources holds the database.
    ch_database = ch_prebuilt_db
        .concat(PHILOSOPHERDATABASE.out.fasta)
        .first()

    emit:
    fasta      = ch_database                                       // channel: [ val(meta), path(fasta) ] - the run database
    fasta_path = ch_database.map { _meta, fasta -> fasta }         // channel: path(fasta) - value channel for broadcasting
}

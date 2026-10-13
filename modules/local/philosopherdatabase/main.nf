include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

/*
 * PHILOSOPHERDATABASE: FragPipe's database update (CmdDatabaseUpdate)
 *
 * FragPipe headless searches database.db-path as given; this update only runs when the user
 * asks for it (--add_philosopher_decoy / --add_philosopher_contaminants), as the FragPipe GUI's "Add decoys" button does.
 * argv: philosopher database --custom <fasta> [--contam --contamprefix] (CmdDatabaseUpdate.java:44-53);
 * the contaminant flags come from ext.args. Philosopher writes its default rev_ decoys.
 *
 * Input:
 *   - meta: Database metadata map
 *   - fasta: Target protein FASTA
 *
 * Output:
 *   - fasta: <prefix>_philosopher.fasta, the FASTA with decoys (and contaminants)
 *   - log: philosopher_database.log, the command line followed by the tool output
 *   - versions: philosopher and FragPipe versions (topic: versions)
 */
process PHILOSOPHERDATABASE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_philosopher.fasta"), emit: fasta
    tuple val(meta), path("philosopher_database.log")   , emit: log
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    export HOME=\$(pwd)

    PHILOSOPHER=${shellQuote(tools.philosopher)}
    [ -x "\$PHILOSOPHER" ] || { echo "ERROR: \$PHILOSOPHER not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }

    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version

    # Philosopher commands run inside a workspace (CmdDatabaseUpdate.usesPhi())
    "\$PHILOSOPHER" workspace --init --nocheck

    CMD=("\$PHILOSOPHER" database --custom ${fasta} ${args})
    printf '%s\\n' "\${CMD[*]}" > philosopher_database.log
    "\${CMD[@]}" 2>&1 | tee -a philosopher_database.log
    mv *-${fasta}.fas ${prefix}_philosopher.fasta

    "\$PHILOSOPHER" workspace --clean --nocheck
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    touch ${prefix}_philosopher.fasta
    echo "philosopher database --custom ${fasta} ${args}" > philosopher_database.log
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version
    """
}

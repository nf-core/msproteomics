include { fragpipeTools; fpopCoadaptrStaging; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FPOPCOADAPTR {
    tag "${meta.id}"
    label 'process_low'

    input:
    // FPOP coADAPTr converter (CmdFPOPcoadaptr.java), run once over the whole FragPipe work dir.
    // inputs: what FragPipe's Files.walk finds in the work dir, by layout
    //   'tree'   - IonQuant results directories, copied whole (every psm.tsv inside is edited)
    //   'groups' - Philosopher per-group result directories ('workdir-root' = blank group);
    //              only their psm.tsv is copied
    //   'diann'  - DIA-NN dia-quant-output directories; only report.tsv is copied and
    //              FPOP_DiannReport_Editor runs instead of FPOP_PSMwriter (CmdFPOPcoadaptr.java:429-433)
    // fpop_mods: comma-joined %.4f FPOP masses (FpopCoadaptrPanel.getFpopMods())
    tuple val(meta), path(inputs, stageAs: 'input/?/*'), val(layout), val(fpop_mods)

    output:
    // Files only: publishing skips files inside a published directory, and
    // these must be published one by one at the replaced task's paths (conf/modules.config).
    tuple val(meta), path("${prefix}/**", type: 'file')                          , emit: tables
    tuple val(meta), path("${prefix}/**/psm.tsv")                                , emit: psms, optional: true
    tuple val(meta), path("${prefix}/**/report.tsv")                             , emit: diann_report, optional: true
    tuple val(meta), path("fpop_coadaptr.log")                                   , emit: log
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version") , emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // FPOP_PSMwriter reads exactly two positional arguments; task.ext.args has no place to go.
    prefix = task.ext.prefix ?: "${meta.id}"
    def editor = layout == 'diann' ? 'FPOP_DiannReport_Editor' : 'FPOP_PSMwriter'
    def tools = fragpipeTools(task)
    // CmdFPOPcoadaptr.java:29-58: java -cp <FragPipe>/lib/* <editor class> <work dir> <FPOP masses>
    def argv = ['java', '-Duser.home={TASK_DIR}', '-cp', "${tools.fragpipe_lib}/*", "org.nesvilab.fragpipe.tools.fpop.${editor}", "{TASK_DIR}/${prefix}", fpop_mods]
    """
    export HOME=\$(pwd)
    echo "${tools.fragpipe_version}" > .fragpipe_version

${fpopCoadaptrStaging(layout, inputs, prefix)}

    CMD=(${shellArgs(argv)})
    printf '%q ' "\${CMD[@]}" > fpop_coadaptr.log
    printf '\\n' >> fpop_coadaptr.log
    "\${CMD[@]}" 2>&1 | tee -a fpop_coadaptr.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def editor = layout == 'diann' ? 'FPOP_DiannReport_Editor' : 'FPOP_PSMwriter'
    def tools = fragpipeTools(task)
    """
    echo "${tools.fragpipe_version}" > .fragpipe_version

    # Same staging as the real run; the psm.tsv files are left unedited.
${fpopCoadaptrStaging(layout, inputs, prefix)}
    printf '%q ' ${shellArgs(['java', '-cp', 'lib/*', "org.nesvilab.fragpipe.tools.fpop.${editor}", prefix, fpop_mods])} > fpop_coadaptr.log
    printf '\\n' >> fpop_coadaptr.log
    """
}

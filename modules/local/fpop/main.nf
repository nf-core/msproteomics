include { fragpipeTools; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FPOP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // FPOP Quant (CmdFpopQuant.java), run once over the whole FragPipe work dir.
    // modpep_tsv:       IonQuant combined_modified_peptide.tsv (label-free mode), [] in TMT mode
    // tmt_report:       TMT-Integrator tmt-report/ directory (TMT mode), [] otherwise
    // tmt_report_unmod: second TMT-Integrator run, tmt-report-unmod/ (TMT mode), [] otherwise
    // config_cli:       FragPipe_FPOP_Analysis.py arguments from parse_fragpipe_workflow.py (CmdFpopQuant.java:90-99),
    //                   naming the inputs as staged here (combined_modified_peptide.tsv, tmt-report/, tmt-report-unmod/)
    tuple val(meta), path(modpep_tsv, stageAs: 'combined_modified_peptide.tsv'), path(tmt_report, stageAs: 'tmt-report'), path(tmt_report_unmod, stageAs: 'tmt-report-unmod'), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                                           , emit: results_dir
    tuple val(meta), path("${prefix}/FPOP_peptides.tsv")                         , emit: peptides
    tuple val(meta), path("${prefix}/FPOP_sites.tsv")                            , emit: sites
    tuple val("${task.process}"), val('fpop'), eval("cat .fpop_version")         , emit: versions_fpop, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version") , emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // FragPipe_FPOP_Analysis.py takes only positional arguments (no flags), so task.ext.args
    // has nowhere to go: an extra word would become the unmodified-peptide argument.
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    if (config_cli[0] == 'combined_modified_peptide.tsv' && !modpep_tsv) {
        error("label-free FPOP needs IonQuant's combined_modified_peptide.tsv")
    }
    // TMT mode reads the TMT-Integrator tables staged as tmt-report/ and tmt-report-unmod/, so the
    // script's outputs (two levels above the multi-site table) land in the task directory
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    echo "${tools.fragpipe_version}" > .fpop_version
    echo "${tools.fragpipe_version}" > .fragpipe_version

    CMD=(${shellArgs(['python3', tools.fpop_script] + config_cli)})
    printf '%q ' "\${CMD[@]}" > ${prefix}/fpop.log
    printf '\\n' >> ${prefix}/fpop.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/fpop.log

    mv FPOP_peptides.tsv FPOP_sites.tsv ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix}
    # Same argument list as the real run; the analysis itself is replaced by empty tables.
    printf '%q ' ${shellArgs(['FragPipe_FPOP_Analysis.py'] + config_cli)} > ${prefix}/fpop.log
    printf '\\n' >> ${prefix}/fpop.log
    touch ${prefix}/FPOP_peptides.tsv ${prefix}/FPOP_sites.tsv
    echo "${tools.fragpipe_version}" > .fpop_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}

include { fragpipeTools; metaproteomicsArgv; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process METAPROTEOMICS {
    tag "${meta.id}"
    label 'process_high'

    input:
    // FP-Meta DbOptimizer (CmdMetaproteomics.java), run once over the whole FragPipe work dir.
    // project_files: every MSFragger <run>.pin and Percolator <run>_percolator_{target,decoy}_psms.tsv,
    //                staged flat like FragPipe's work dir (FP-Meta walks it 2 levels deep)
    // config_cli:    DbOptimizer values (qvalue, deltaHyperscore, minPeptCntPerProt, minUniqPeptCntPerProt,
    //                minUniqPeptCnt, hostName, iterations, cmdLineOpts) from parse_fragpipe_workflow.py
    tuple val(meta), path(project_files, stageAs: 'project/*'), val(config_cli)
    path(fasta)            // database FASTA file (shared resource)
    path(taxon_name_file)  // NCBI names.dmp; [] = FragPipe's tools/metaproteomics/names.dmp
    path(taxon_node_file)  // NCBI nodes.dmp; [] = FragPipe's tools/metaproteomics/nodes.dmp
    val decoy_tag          // decoy protein prefix from workflow database.decoy-tag (FragpipeRun.java:1442)

    output:
    tuple val(meta), path("${prefix}")                                                       , emit: results_dir
    tuple val("${task.process}"), val('metaproteomics'), eval("cat .metaproteomics_version") , emit: versions_metaproteomics, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version")             , emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools = fragpipeTools(task)
    def inputs = (project_files instanceof List ? project_files : [project_files]).collect { f -> f.name.tokenize('/').last() }
    // CmdMetaproteomics.java:39-42: the taxonomy dumps ship in FragPipe's tools/metaproteomics/
    def names_dmp = taxon_name_file ? "{TASK_DIR}/${taxon_name_file}" : tools.taxon_names
    def nodes_dmp = taxon_node_file ? "{TASK_DIR}/${taxon_node_file}" : tools.taxon_nodes
    // FP-Meta walks the project dir without FOLLOW_LINKS, so the staged links are copied; FragPipe passes its
    // work dir as both --projectDir and --outdir (CmdMetaproteomics.java:75-80), so FP-Meta's outputs are
    // that dir without the inputs. argv of CmdMetaproteomics.java:70-111, cmd-line-opts one element.
    def argv = ['java', '-Duser.home={TASK_DIR}', '-Djava.awt.headless=true', "-Xmx${ram}G", '-jar', tools.fpmeta_jar] +
        metaproteomicsArgv(config_cli, '{TASK_DIR}/project_resolved', "{TASK_DIR}/${fasta}", decoy_tag, names_dmp, nodes_dmp)
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    echo "${tools.fpmeta_version}" > .metaproteomics_version
    echo "${tools.fragpipe_version}" > .fragpipe_version

    cp -rL project project_resolved

    CMD=(${shellArgs(argv)} ${args})
    printf '%q ' "\${CMD[@]}" > ${prefix}/metaproteomics.log
    printf '\\n' >> ${prefix}/metaproteomics.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/metaproteomics.log

    rm ${inputs.collect { n -> "'project_resolved/${n}'" }.join(' ')}
    mv project_resolved/* ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def inputs = (project_files instanceof List ? project_files : [project_files]).collect { f -> f.name.tokenize('/').last() }
    """
    mkdir -p ${prefix}
    printf '%q ' ${shellArgs(['DbOptimizer', '--decoyTag', decoy_tag, '--qvalue', config_cli.qvalue, '--hostName', config_cli.hostName, '--iterations', config_cli.iterations])} > ${prefix}/metaproteomics.log
    printf '\\n' >> ${prefix}/metaproteomics.log
    printf '%s\\n' ${shellArgs(inputs.sort())} >> ${prefix}/metaproteomics.log
    touch ${prefix}/refined_db.fasta ${prefix}/combined_PSM_species.tsv
    echo "${tools.fpmeta_version}" > .metaproteomics_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}

/*
 * MERGESPLITSEARCH: Merge MSFragger split-database search results per sample.
 *
 * For each sample, collects pepXML, PIN, and score histogram files from all
 * database chunks, sums histograms, generates expect functions via MSFragger,
 * and merges search results by re-ranking hits per spectrum.
 *
 * Faithfully ports the merge logic from FragPipe's msfragger_pep_split.py.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - chunk_dirs: the MSFRAGGER --partial result directories of the sample, in chunk order
 *     (pepXML, PIN and score histogram per chunk), staged as chunk_1..chunk_N
 *   - params_text: text of the file's search params (output_report_topN, output_max_expect,
 *     output_file_extension; msfragger_pep_split.py:76-80)
 *   - fasta: the full FASTA; its path in the task is written into the merged pepXML header
 *
 * Output:
 *   - pepxml: Merged pepXML file
 *   - pin: Merged PIN file
 *   - versions: Software versions
 */
include { fragpipeJavaMemGb; fragpipeTools; requireFragpipeLicenseAgreement; writeTextFile; shellArgs; msfraggerSplitMergeSettings; fragpipeJarLines } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process MERGESPLITSEARCH {
    tag "${meta.id}"
    label 'process_medium'

    // No container directive — inherits process.container (the FragPipe image)
    // because this module uses MSFragger JAR for --generate_expect_functions

    input:
    tuple val(meta), path(chunk_dirs, stageAs: 'chunk_*'), val(params_text)
    // fasta: the full FASTA. msfragger_pep_split.py writes the params' database_name, the FASTA path the tools
    //   run with, into the merged header (:92-94, :302-303); here its path in the task, as MSFRAGGER's database_name
    path fasta

    output:
    tuple val(meta), path("${prefix}/*.pepXML"), emit: pepxml
    tuple val(meta), path("${prefix}/*.pin")   , emit: pin
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def agreement = requireFragpipeLicenseAgreement(task)
    def mem = fragpipeJavaMemGb(task)
    def num_chunks = (chunk_dirs instanceof List ? chunk_dirs : [chunk_dirs]).size()
    def merge = msfraggerSplitMergeSettings(params_text)
    // msfragger_cmd of msfragger_pep_split.py:73 (jvm_cmd + jar); generate_expect_cmd adds --license <file> (:114).
    // MSFragger finds ext/thermo and ext/bruker next to its jar; FragPipe passes no -Dlibs.* to it.
    def msfragger_cmd = "java -Duser.home={TASK_DIR} -jar -Dfile.encoding=UTF-8 -Xmx${mem}G ${tools.msfragger_jar} --license ${tools.license}"
    """
    ${writeTextFile('I_AGREE_FRAGPIPE_LICENSE_AGREEMENT', agreement)}
    ${fragpipeJarLines(tools, ['msfragger'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version

    merge_split_search.py \\
        --sample_name ${prefix} \\
        --chunk_dirs ${(1..num_chunks).collect { i -> "chunk_${i}" }.join(',')} \\
        --num_chunks ${num_chunks} \\
        --msfragger_cmd ${shellArgs([msfragger_cmd])} \\
        --outdir . \\
        --output_report_topN ${merge.output_report_topN} \\
        --output_max_expect ${merge.output_max_expect} \\
        --output_file_extension ${merge.output_file_extension} \\
        --fasta_path ${shellArgs(["{TASK_DIR}/${fasta}"])} ${args}

    # Merged next to the LC-MS file as FragPipe writes them (summary_xml, msfragger_pep_split.py:397, :462),
    # then into ${prefix}/ like MSFRAGGER's results
    mkdir -p ${prefix}
    mv ${prefix}*.${merge.output_file_extension} ${prefix}.pin ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    mkdir -p "${prefix}"
    touch "${prefix}/${prefix}.pepXML"
    touch "${prefix}/${prefix}.pin"
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    echo "${tools.fragpipe_version}" > .msfragger_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}

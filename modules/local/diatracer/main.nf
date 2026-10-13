include { fragpipeTools; requireFragpipeLicenseAgreement; writeTextFile; fragpipeJarLines; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process DIATRACER {
    tag "${meta.id}"
    label 'process_high'

    input:
    // diaTracer processes Bruker .d files for DIA analysis
    // d_file: Bruker .d directory
    // config_cli: diaTracer parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(d_file), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/*_diatracer.mzML")         , emit: mzml
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('diatracer'), eval("cat .diatracer_version"), emit: versions_diatracer, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools = fragpipeTools(task)
    def agreement = requireFragpipeLicenseAgreement(task)
    // Argv as CmdDiaTracer.java:84-117 (headless run): the Bruker native libraries next to the MSFragger jar
    // (CmdMsfragger.searchExtLibsBruker, CmdMsfragger.java:264-273), the license FragPipe finds in tools/..
    // (FragpipeLocations.java:340-353, passed by path only); config_cli holds --writeInter ... --RFMax,
    // then ext.args: tool arguments, not JVM options
    // CmdDiaTracer.java:138: diaTracer writes <the .d path with ".d" replaced by "_diatracer.mzML"> next to the .d it is
    // given, so the staged path (not its resolved target in the data) is passed and the mzML lands in the task dir
    def mzml_name = d_file.name.replace('.d', '_diatracer.mzML')
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    ${writeTextFile('I_AGREE_FRAGPIPE_LICENSE_AGREEMENT', agreement)}
    ${fragpipeJarLines(tools, ['diatracer', 'msfragger'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version

    CMD=(java -Duser.home=\$PWD -Djava.awt.headless=true -Xmx${ram}G ${shellArgs(["-Dlibs.bruker.dir=${tools.bruker_dir}"])} -jar ${shellArgs([tools.diatracer_jar])} --license "${tools.license}" --dFilePath "\$(pwd)/${d_file}" --workDir "\$(pwd)/${prefix}" --threadNum ${task.cpus} ${config_cli} ${args})
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/diatracer.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/diatracer.log
    mv "${mzml_name}" ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def mzml_name = d_file.name.replace('.d', '_diatracer.mzML')
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix}
    touch "${prefix}/${mzml_name}"
    touch ${prefix}/diatracer.log
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    echo "${tools.diatracer_version}" > .diatracer_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}

process SKYLINE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // Skyline creates Skyline documents from FragPipe results
    // results_dir: FragPipe output directory containing psm.tsv, library files
    // mzml_files: mass spec files for import
    // speclib: spectral library file
    // config_cli: Skyline parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(results_dir), path(mzml_files, stageAs: 'spectra/*'), path(speclib), val(config_cli)
    val(skyline_path)   // Path to Skyline executable (shared resource)

    output:
    tuple val(meta), path("${prefix}")                              , emit: results_dir
    tuple val(meta), path("${prefix}/skyline_files/fragpipe.sky")   , emit: skyline_document, optional: true
    tuple val(meta), path("${prefix}/*.csv")                        , emit: reports, optional: true
    tuple val("${task.process}"), val('skyline'), eval("cat .skyline_version"), emit: versions_skyline, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // Skyline module (via FragPipe Java wrapper) CLI parameters:
    // skylinePath: path to Skyline executable
    // workDir: output directory
    // skylineVersion: Skyline version string
    // modsMode: modification handling mode (0-3)
    // precursorTolerance: precursor mass tolerance in ppm
    // fragmentTolerance: fragment mass tolerance in ppm
    // rtTolerance: retention time tolerance in minutes
    // libraryProductIons: number of product ions from library
    // runSkylineQuant: run Skyline quantification

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}/skyline_files

    TOOLS_DIR="${tools_dir}"

    # JAR discovery: fragpipe JAR + batmass-io + fragpipe-lib
    FP_JAR=\$(find "\$TOOLS_DIR" -path '*/lib/fragpipe*.jar' -type f 2>/dev/null | head -1 || true)
    BATMASS_JAR=\$(find "\$TOOLS_DIR" -name 'batmass-io*.jar' -type f 2>/dev/null | head -1 || true)
    FP_LIB="${tools_dir}/../lib"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .skyline_version

    echo "\$_fp_ver" > .fragpipe_version

    # Run Skyline via FragPipe Java wrapper
    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \$JAVA_OPTS -cp \\"\$FP_JAR:\$BATMASS_JAR:\$FP_LIB/*\\" org.nesvilab.fragpipe.tools.skyline.Skyline ${skyline_path} ${prefix} auto ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/skyline.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/skyline.log

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/skyline_files
    touch ${prefix}/skyline_files/fragpipe.sky
    touch ${prefix}/skyline.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .skyline_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}

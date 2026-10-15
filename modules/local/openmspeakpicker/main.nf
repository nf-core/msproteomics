// Derived from nf-core/quantms OPENMSPEAKPICKER module
// Changes: options via ext.args, output keeps the input file name, added stub/when clause, updated container
process OPENMSPEAKPICKER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84ad5437e080c0b3ec0ce315a54c0bb62950277a4496f4d1be911f82078dd2a1/data' :
        'community.wave.seqera.io/library/openms:3.6.0--88df48c0a426b7f7' }"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("picked/${mzml_file.name}"), emit: mzmls_picked
    tuple val(meta), path("${prefix}_pp.log")         , emit: log
    tuple val("${task.process}"), val('openms'), eval("PeakPickerHiRes 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p picked
    PeakPickerHiRes \\
        -in ${mzml_file} \\
        -out picked/${mzml_file.name} \\
        -threads $task.cpus \\
        $args \\
        2>&1 | tee ${prefix}_pp.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p picked
    echo "<mzML/>" > picked/${mzml_file.name}
    echo "Stub execution" > ${prefix}_pp.log
    """
}

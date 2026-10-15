// Derived from nf-core/quantms OPENMSPEAKPICKER module
// Changes: options via ext.args, output keeps the input file name, added stub/when clause, updated container
process OPENMSPEAKPICKER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

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

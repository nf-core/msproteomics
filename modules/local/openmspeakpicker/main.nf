// Derived from nf-core/quantms OPENMSPEAKPICKER module
// Changes: options via ext.args, output keeps the input file name, added stub/when clause, updated container
process OPENMSPEAKPICKER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("picked/${mzml_file.name}"), emit: mzmls_picked
    tuple val("${task.process}"), val('openms'), eval("cat .openms_peakpicker_version"), emit: versions_openms, topic: versions
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p picked
    PeakPickerHiRes \\
        -in ${mzml_file} \\
        -out picked/${mzml_file.name} \\
        -threads $task.cpus \\
        $args \\
        2>&1 | tee ${prefix}_pp.log

    PeakPickerHiRes 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1 > .openms_peakpicker_version
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p picked
    echo "<mzML/>" > picked/${mzml_file.name}
    echo "Stub execution" > ${prefix}_pp.log
    echo "stub" > .openms_peakpicker_version
    """
}

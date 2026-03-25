// Derived from nf-core/quantms OPENMSPEAKPICKER module
// Changes: fixed "inmermory" typo, params moved to task.ext, added stub/when clause, updated container
process OPENMS_PEAK_PICKER {
    tag "$meta.mzml_id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/bigbio/openms-tools-thirdparty-sif:2025.04.14' :
        'ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14' }"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*.mzML"), emit: mzmls_picked
    tuple val("${task.process}"), val('openms'), eval("cat .openms_peakpicker_version"), emit: versions_openms, topic: versions
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    def in_mem = task.ext.peakpicking_inmemory ? "inmemory" : "lowmemory"
    def lvls = task.ext.peakpicking_ms_levels ? "-algorithm:ms_levels ${task.ext.peakpicking_ms_levels}" : ""
    def pp_debug = task.ext.pp_debug ?: 0

    """
    PeakPickerHiRes \\
        -in ${mzml_file} \\
        -out ${prefix}.mzML \\
        -threads $task.cpus \\
        -debug $pp_debug \\
        -processOption ${in_mem} \\
        ${lvls} \\
        $args \\
        2>&1 | tee ${prefix}_pp.log

    PeakPickerHiRes 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1 > .openms_peakpicker_version
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    echo "<mzML/>" > ${prefix}.mzML
    echo "Stub execution" > ${prefix}_pp.log
    echo "stub" > .openms_peakpicker_version
    """
}

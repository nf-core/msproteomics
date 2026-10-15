/*
 * MZMLINDEXING: Rewrite an mzML file as indexed mzML (OpenMS FileConverter)
 *
 * Input:
 *   - meta: Sample metadata map
 *   - mzmlfile: mzML file
 *
 * Output:
 *   - mzmls_indexed: out/<name>.mzML
 *   - versions: versions.yml (topic: versions)
 *   - log: FileConverter output
 */
process MZMLINDEXING {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(mzmlfile)

    output:
    tuple val(meta), path("out/*.mzML"), emit: mzmls_indexed
    path "versions.yml", emit: versions, topic: versions
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p out
    FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML $args 2>&1 | tee ${mzmlfile.baseName}_mzmlindexing.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FileConverter: \$(FileConverter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    touch out/${mzmlfile.baseName}.mzML
    echo "Stub execution" > ${mzmlfile.baseName}_mzmlindexing.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FileConverter: \$(FileConverter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """
}

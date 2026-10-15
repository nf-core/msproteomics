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
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84ad5437e080c0b3ec0ce315a54c0bb62950277a4496f4d1be911f82078dd2a1/data' :
        'community.wave.seqera.io/library/openms:3.6.0--88df48c0a426b7f7' }"

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

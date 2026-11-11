include { textFileCommand } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/write_text_file'

/*
 * WRITETEXTFILE: Write a text file whose body was built outside the tool modules (Groovy helpers or the
 * FragPipe's plan), so the tool modules receive it as a path and only run their tool.
 *
 * Input:
 *   - meta: Metadata map carried to the output
 *   - name: File name to write
 *   - text: File body, written byte for byte (a body that does not end with a newline gets none)
 *
 * Output:
 *   - file: The written file
 *   - versions_coreutils: coreutils version (cat)
 */
process WRITETEXTFILE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta), val(name), val(text)

    output:
    tuple val(meta), path("${name}"), emit: file
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed -n '1s/.* //p'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def write = textFileCommand(name, text)
    """
    ${write}
    """

    stub:
    def write = textFileCommand(name, text)
    """
    ${write}
    """
}

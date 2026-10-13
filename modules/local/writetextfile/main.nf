include { textFileCommand } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/write_text_file'

/*
 * WRITETEXTFILE: Write a text file whose body was built outside the tool modules (Groovy helpers or the
 * workflow parser), so the tool modules receive it as a path and only run their tool.
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
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

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

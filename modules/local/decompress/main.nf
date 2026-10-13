
process DECOMPRESS {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' :
        'quay.io/biocontainers/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' }"

    stageInMode {
        if (task.attempt == 1) {
            if (task.executor == 'awsbatch') {
                'symlink'
            } else {
                'link'
            }
        } else if (task.attempt == 2) {
            if (task.executor == 'awsbatch') {
                'copy'
            } else {
                'symlink'
            }
        } else {
            'copy'
        }
    }

    input:
    tuple val(meta), path(compressed_file)

    output:
    tuple val(meta), path("${target_name}"), emit: decompressed_files
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version 2>&1 | head -1 | grep -oE "[0-9]+\\.[0-9]+(\\.[0-9]+)?" | head -1'), emit: versions_gunzip, topic: versions
    tuple val("${task.process}"), val('tar'), eval('tar --version 2>&1 | head -1 | grep -oE "[0-9]+\\.[0-9]+(\\.[0-9]+)?" | head -1'), emit: versions_tar, topic: versions
    tuple val("${task.process}"), val('unzip'), eval('unzip -v 2>&1 | head -1 | grep -oE "[0-9]+\\.[0-9]+" | head -1'), emit: versions_unzip, topic: versions
    path '*.log'                           , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // X.d.tar.gz -> X.d, X.mzML.gz -> X.mzML, X.raw.zip -> X.raw
    target_name = compressed_file.name.replaceAll(/\.(tar\.gz|tar|gz|zip)$/, '')
    """
    case ${compressed_file} in
        *.tar.gz) tar xvzf ${compressed_file} ;;
        *.tar)    tar xvf ${compressed_file} ;;
        *.gz)     gunzip -c ${compressed_file} > ${target_name} ;;
        *.zip)    unzip ${compressed_file} ;;
        *)        echo "Unknown archive type: ${compressed_file}" >&2; exit 1 ;;
    esac 2>&1 | tee ${prefix}_decompression.log

    # A Bruker archive may hold the .d directory under another name.
    if [ ! -e ${target_name} ]; then
        found=\$(find . -maxdepth 1 -name '*.d' -type d ! -name '.')
        if [ "\$(printf '%s\\n' "\$found" | grep -c .)" != 1 ]; then
            echo "Archive ${compressed_file} did not contain ${target_name} (found: \$found)" >&2
            exit 1
        fi
        mv "\$found" ${target_name}
    fi

    ls -l | tee -a ${prefix}_decompression.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    target_name = compressed_file.name.replaceAll(/\.(tar\.gz|tar|gz|zip)$/, '')
    def create_target = target_name.endsWith('.d') ? "mkdir -p ${target_name} && touch ${target_name}/analysis.tdf" : "touch ${target_name}"
    """
    ${create_target}
    echo "Stub execution" > ${prefix}_decompression.log
    """
}

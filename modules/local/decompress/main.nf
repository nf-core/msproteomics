
process DECOMPRESS {
    tag "$meta.mzml_id"
    label 'process_low'
    label 'error_retry'
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
    tuple val(meta), path('*.d'),   emit: decompressed_files
    path 'versions.yml',   emit: versions
    path '*.log',   emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    String prefix = task.ext.prefix ?: "${meta.mzml_id}"
    def target_name = file(compressed_file.baseName).baseName

    """
    function verify_tar {
        exit_code=0
        error=\$(tar df \$1 2>&1) || exit_code=\$?
        if [ \$exit_code -eq 2 ]; then
            echo "\${error}"
            exit 2
        fi

        case \${error} in
            *'No such file'* )
                echo "\${error}" | grep "No such file"
                exit 1
                ;;
            *'Size differs'* )
                echo "\${error}" | grep "Size differs"
                exit 1
                ;;
        esac
    }

    function extract {
        if [ -z "\$1" ]; then
            echo "Usage: extract <path/file_name>.<gz|tar|tar.bz2>"
        else
            if [ -f \$1 ]; then
                case \$1 in
                    *.tar.gz)    tar xvzf \$1 && verify_tar \$1               ;;
                    *.gz)        gunzip \$1                                     ;;
                    *.tar)       tar xvf \$1 && verify_tar \$1                ;;
                    *.zip)       unzip \$1                                      ;;
                    *)           echo "extract: '\$1' - unknown archive method" ;;
                esac
            else
                echo "\$1 - file does not exist"
            fi
        fi
    }

    tar --help 2>&1 | tee -a ${prefix}_decompression.log
    gunzip --help 2>&1 | tee -a ${prefix}_decompression.log
    (unzip --help 2>&1 || zip --help 2>&1) | tee -a ${prefix}_decompression.log
    echo "Unpacking..." | tee -a ${prefix}_decompression.log

    extract ${compressed_file} 2>&1 | tee -a ${prefix}_decompression.log
    [ -d ${target_name}.d ] && \\
        echo "Found ${target_name}.d" || \\
        mv *.d ${target_name}.d

    ls -l | tee -a ${prefix}_decompression.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(gunzip --help 2>&1 | head -1 | grep -oE "\\d+\\.\\d+(\\.\\d+)?")
        tar: \$(tar --help 2>&1 | head -1 | grep -oE "\\d+\\.\\d+(\\.\\d+)?")
        unzip: \$((unzip --help 2>&1 || zip --help 2>&1) | head -2 | tail -1 | grep -oE "\\d+\\.\\d+")
    END_VERSIONS
    """

    stub:
    String prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    mkdir -p ${prefix}.d
    echo "Stub execution" > ${prefix}_decompression.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: stub
        tar: stub
        unzip: stub
    END_VERSIONS
    """
}

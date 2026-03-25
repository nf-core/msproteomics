process TMTINTEGRATOR {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(psm_dirs, stageAs: 'input_results/*')
    path(annotation_files)
    path(config_file)

    output:
    tuple val(meta), path("${prefix}")                     , emit: results_dir
    tuple val(meta), path("${prefix}/abundance_*.tsv")     , emit: abundance, optional: true
    tuple val(meta), path("${prefix}/ratio_*.tsv")         , emit: ratio, optional: true
    tuple val(meta), path("${prefix}/*_Proteins.txt")      , emit: proteins_txt, optional: true
    tuple val(meta), path("${prefix}/*_Peptides.txt")      , emit: peptides_txt, optional: true
    tuple val("${task.process}"), val('tmtintegrator'), eval("cat .tmtintegrator_version"), emit: versions_tmtintegrator, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // CRITICAL: No 'def' for prefix — must be visible in output block (Nextflow 25.x scoping)
    prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    TMTINTEGRATOR_JAR=\$(find "\$TOOLS_DIR" -name 'TMT-Integrator*.jar' -type f 2>/dev/null | head -1 || true)
    _ver=\$(java -jar "\$TMTINTEGRATOR_JAR" 2>&1 | grep -oP 'TMT Integrator \\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .tmtintegrator_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Copy config file (if provided) and append runtime-specific values
    if [ -s "${config_file}" ]; then
        cp ${config_file} ${prefix}/tmt-integrator-conf.yml
    else
        # No config provided — create minimal YAML with runtime params only;
        # TMTIntegrator will use its built-in defaults for algorithm parameters.
        touch ${prefix}/tmt-integrator-conf.yml
    fi
    cat >> ${prefix}/tmt-integrator-conf.yml <<EOF
  # Runtime values (appended by module)
  path: \$TMTINTEGRATOR_JAR
  memory: ${mem}
  output: \$(pwd)/${prefix}
EOF

    # TMTIntegrator determines plexes from directory structure:
    #   each psm.tsv path → parent directory name → experiment → plex
    # All LCMS files from the same plex must be in ONE directory.
    #
    # Annotation format detection:
    #   6+ columns (tab-separated header starts with "plex") → TMTIntegrator format (multi-plex)
    #   2 columns (space/tab-separated, no header) → IonQuant format (single-plex)
    #
    # Single-plex: combine all per-sample psm.tsv into one experiment directory
    # Multi-plex: pass per-sample directories directly (each sample = one plex)

    MULTI_PLEX=false
    if [[ -f "${annotation_files}" && -s "${annotation_files}" ]]; then
        FIRST_LINE=\$(head -1 "${annotation_files}")
        if echo "\$FIRST_LINE" | grep -qP '^plex\\t'; then
            MULTI_PLEX=true
            cp "${annotation_files}" experiment_annotation.tsv
        fi
    fi

    if [[ "\$MULTI_PLEX" == "false" ]]; then
        # Single-plex mode: combine all per-sample psm.tsv and protein.tsv into
        # ONE experiment directory (CmdPhilosopherReport.java combined output pattern).
        #
        # Input directories are staged under input_results/ via stageAs to avoid
        # name collision with the output prefix (e.g., IONQUANT_ISOBARIC outputs
        # a directory named "${prefix}" which would clash with our combined dir).
        mkdir -p "${prefix}"

        # Find source psm.tsv at any depth within the staged input_results/ directory.
        # Concatenate psm.tsv: header from first file, data from all files
        FIRST=true
        for psm_file in \$(find -L input_results -name "psm.tsv" | sort); do
            if \$FIRST; then
                cat "\$psm_file" > "${prefix}/psm.tsv"
                FIRST=false
            else
                tail -n +2 "\$psm_file" >> "${prefix}/psm.tsv"
            fi
        done

        if \$FIRST; then
            echo "ERROR: No psm.tsv files found in staged input directories."
            echo "This usually means the channel routing did not pass the correct directories to TMTIntegrator."
            exit 1
        fi

        # Validate that psm.tsv contains TMT intensity columns (columns starting with "Intensity ")
        if ! head -1 "${prefix}/psm.tsv" | grep -q "Intensity "; then
            echo "ERROR: psm.tsv lacks TMT intensity columns (no 'Intensity ' columns found)."
            echo "Header: \$(head -1 "${prefix}/psm.tsv")"
            echo "Column count: \$(head -1 "${prefix}/psm.tsv" | awk -F'\\t' '{print NF}')"
            echo "This usually means IONQUANT_ISOBARIC or PHILOSOPHER_LABELQUANT output was not correctly passed to TMTIntegrator."
            exit 1
        fi

        # Combine protein.tsv: header from first, union of all proteins (deduplicated)
        FIRST=true
        for protein_file in \$(find -L input_results -name "protein.tsv" | sort); do
            if \$FIRST; then
                cat "\$protein_file" > "${prefix}/protein.tsv"
                FIRST=false
            else
                tail -n +2 "\$protein_file" >> "${prefix}/protein.tsv"
            fi
        done
        if [[ -f "${prefix}/protein.tsv" ]]; then
            head -1 "${prefix}/protein.tsv" > "${prefix}/protein_dedup.tsv"
            tail -n +2 "${prefix}/protein.tsv" | sort -t\$'\\t' -k1,1 -u >> "${prefix}/protein_dedup.tsv"
            mv "${prefix}/protein_dedup.tsv" "${prefix}/protein.tsv"
        fi

        # Generate experiment_annotation.tsv from IonQuant annotation (2-col → 6-col)
        if [[ -f "${annotation_files}" && -s "${annotation_files}" ]]; then
            printf 'plex\\tchannel\\tsample\\tsample_name\\tcondition\\treplicate\\n' > experiment_annotation.tsv
            while IFS=\$' \\t' read -r channel sample_name; do
                [[ -z "\$channel" ]] && continue
                condition=\$(echo "\$sample_name" | sed 's/_[0-9]*\$//')
                replicate=\$(echo "\$sample_name" | grep -oP '_\\K[0-9]+\$')
                if [[ -z "\$replicate" ]]; then
                    echo "ERROR: Cannot extract replicate number from sample_name '\$sample_name' (expected trailing _N)" >&2
                    exit 1
                fi
                printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${prefix}" "\$channel" "\$sample_name" "\$condition" "\$condition" "\$replicate" >> experiment_annotation.tsv
            done < "${annotation_files}"
        fi
    fi

    # Collect psm.tsv paths for TMTIntegrator
    if [[ "\$MULTI_PLEX" == "false" ]]; then
        # Single-plex: use only the combined psm.tsv we just created
        PSM_PATHS="./${prefix}/psm.tsv"
    else
        # Multi-plex: find all per-sample psm.tsv within the staged input directory
        PSM_PATHS=\$(find -L input_results -name "psm.tsv" | sort | tr '\\n' ' ')
    fi

    export JAVA_OPTS="-Xmx${mem}G"
    java \${JAVA_OPTS:-} -jar "\$TMTINTEGRATOR_JAR" \\
        ${prefix}/tmt-integrator-conf.yml \\
        \$PSM_PATHS ${args} \\
        2>&1 | tee ${prefix}/tmt-integrator.log

    """

    stub:
    // CRITICAL: No 'def' for prefix — must be visible in output block (Nextflow 25.x scoping)
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/abundance_gene_None.tsv
    touch ${prefix}/abundance_protein_None.tsv
    touch ${prefix}/ratio_gene_None.tsv
    touch ${prefix}/ratio_protein_None.tsv
    touch ${prefix}/tmt-integrator.log
    touch ${prefix}/tmt-integrator-conf.yml
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .tmtintegrator_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}

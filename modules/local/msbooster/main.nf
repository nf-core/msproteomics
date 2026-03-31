process MSBOOSTER {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pin_files), path(mzml_files)
    path params_file      // Full MSBooster params file (key = value). Pass [] when not using.
    path fragger_params   // MSFragger params file for modification definitions. Pass [] when not using.
                          // FragPipe writes "fragger = /path/to/fragger.params" (CmdMSBooster.java:171)
                          // so MSBooster can read mod definitions for rescoring.
    val has_ion_mobility  // true if input data has ion mobility (e.g., Bruker timsTOF .d files).
                          // Used to auto-detect useIM when not set in params file
                          // (matching FragPipe CmdMSBooster.java:152-170 auto-detection).

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val(meta), path("${prefix}/*.pin")        , emit: pin_edited
    tuple val("${task.process}"), val('msbooster'), eval("cat .msbooster_version"), emit: versions_msbooster, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def pin_list = pin_files instanceof List ? pin_files.join(' ') : pin_files.toString()
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    MSBOOSTER_JAR=\$(find "\$TOOLS_DIR" -name 'MSBooster*.jar' -type f 2>/dev/null | head -1 || true)
    BATMASS_JAR=\$(find "\$TOOLS_DIR" -name 'batmass-io*.jar' -type f 2>/dev/null | head -1 || true)
    # DIA-NN binary: used by MSBooster for spectral/RT prediction
    DIANN_BIN=\$(find "\$TOOLS_DIR" -path "*/diann/*/linux/diann*" -type f 2>/dev/null | sort | tail -1)
    if [ -z "\$DIANN_BIN" ]; then
        # Fall back to PATH (e.g., commercial container with diann in /usr/local/bin)
        DIANN_BIN="diann"
    fi
    # unimod.obo: required by MSBooster for modification definitions
    UNIMOD_OBO="\$TOOLS_DIR/unimod.obo"
    if [ ! -f "\$UNIMOD_OBO" ]; then
        UNIMOD_OBO="/usr/local/share/unimod.obo"
    fi
    _ver=\$(java -jar "\$MSBOOSTER_JAR" --version 2>&1 | grep -oP 'v\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .msbooster_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    if [ -s "${params_file}" ]; then
        # Full params file mode: matches FragPipe calling convention
        # Start with user params, then add/override runtime params
        cp ${params_file} ${prefix}/msbooster_params.txt
        sed -i '/^numThreads/d' ${prefix}/msbooster_params.txt
        sed -i '/^DiaNN/d' ${prefix}/msbooster_params.txt
        sed -i '/^unimodObo/d' ${prefix}/msbooster_params.txt
        sed -i '/^fragger /d' ${prefix}/msbooster_params.txt
        sed -i '/^mzmlDirectory/d' ${prefix}/msbooster_params.txt
        sed -i '/^pinPepXMLDirectory/d' ${prefix}/msbooster_params.txt
        echo "numThreads = ${task.cpus}" >> ${prefix}/msbooster_params.txt
        echo "DiaNN = \$DIANN_BIN" >> ${prefix}/msbooster_params.txt
        echo "unimodObo = \$UNIMOD_OBO" >> ${prefix}/msbooster_params.txt
        # fragger param: MSBooster reads modification definitions from fragger.params
        # (CmdMSBooster.java:171). Only write if fragger_params file is provided.
        if [ -s "${fragger_params}" ]; then
            echo "fragger = \$(pwd)/${fragger_params}" >> ${prefix}/msbooster_params.txt
        fi
        echo "mzmlDirectory = \$(pwd)" >> ${prefix}/msbooster_params.txt
        # pinPepXMLDirectory: FragPipe writes relative paths (e.g. sample1/sample1.pin),
        # but MSBooster handles both absolute and relative paths equivalently.
        # Nextflow stages files in the working directory, so flat filenames work correctly.
        echo "pinPepXMLDirectory = ${pin_list}" >> ${prefix}/msbooster_params.txt
    else
        # CLI mode (backward compatible): minimal params file + ext.args
        # Defaults match FragPipe behavior (CmdMSBooster.java:164,167,172)
        cat > ${prefix}/msbooster_params.txt <<PARAMS
useDetect = false
renamePin = 1
deletePreds = false
useRT = true
useSpectra = true
useIM = false
numThreads = ${task.cpus}
DiaNN = \$DIANN_BIN
unimodObo = \$UNIMOD_OBO
mzmlDirectory = \$(pwd)
pinPepXMLDirectory = ${pin_list}
PARAMS
        # fragger param: MSBooster reads modification definitions from fragger.params
        if [ -s "${fragger_params}" ]; then
            echo "fragger = \$(pwd)/${fragger_params}" >> ${prefix}/msbooster_params.txt
        fi
    fi

    # Auto-detect useIM for ion mobility data (matching FragPipe CmdMSBooster.java:152-170).
    # FragPipe enables useIM only when BOTH conditions are true:
    #   1. Input data is timsTOF (hasTimsTof) — detected from .d file extensions
    #   2. User enabled IM prediction (predictIm) — from workflow's msbooster.predict-im
    # If useIM is already in the params file (from workflow's msbooster.predict-im),
    # that explicit value is respected. Otherwise, auto-detect from data type.
    if ! grep -q '^useIM' ${prefix}/msbooster_params.txt; then
        if [ "${has_ion_mobility}" = "true" ]; then
            echo "useIM = true" >> ${prefix}/msbooster_params.txt
        else
            echo "useIM = false" >> ${prefix}/msbooster_params.txt
        fi
    fi

    export JAVA_OPTS="-Xmx${mem}G"
    java \$JAVA_OPTS -cp "\$MSBOOSTER_JAR:\$BATMASS_JAR" mainsteps.MainClass \\
        --paramsList ${prefix}/msbooster_params.txt \\
        ${args} \\
        2>&1 | tee ${prefix}/${prefix}_msbooster.log
    MSBOOSTER_EXIT=\${PIPESTATUS[0]}

    shopt -s nullglob
    for f in *_edited.pin; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob

    # MSBooster appends '_edited' suffix to output PIN files; normalize to standard naming.
    # FragPipe passes _edited.pin to Percolator when MSBooster is enabled (CmdPercolator.java:207-211).
    # Here we strip the suffix so downstream modules (Percolator) receive consistently named .pin files
    # regardless of whether MSBooster was run. The fragpipe_search subworkflow handles routing.
    for f in ${prefix}/*_edited.pin; do
        [ -f "\$f" ] && mv "\$f" "\${f%_edited.pin}.pin"
    done

    # Check if output file exists
    EDITED_PIN_COUNT=\$(find ${prefix} -name "*.pin" -type f | wc -l)
    if [ "\$EDITED_PIN_COUNT" -eq 0 ]; then
        echo "ERROR: MSBooster failed to create edited pin file (exit code: \$MSBOOSTER_EXIT)" >&2
        exit 1
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def pin_input = pin_files instanceof List ? pin_files.join(' ') : pin_files.toString()
    """
    mkdir -p ${prefix}
    # Create output PIN for each input PIN (supports aggregate mode with multiple PINs)
    for f in ${pin_input}; do
        stem=\$(basename "\$f" .pin)
        touch ${prefix}/\${stem}.pin
    done
    touch ${prefix}/${prefix}_msbooster.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .msbooster_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}

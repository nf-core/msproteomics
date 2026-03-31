process IONQUANT {
    tag "${meta.id}"
    label 'process_high'

    input:
    // psm_dirs: PHILOSOPHER_FILTER output directories, each named by sample and containing psm.tsv
    // spec_files: mzML files staged in spectra/ subdirectory
    tuple val(meta), path(psm_dirs), path(spec_files, stageAs: 'spectra/*')
    path(annotation_file)
    val(config_cli)  // Config as CLI string (e.g., "--mbr 1 --maxlfq 1 --perform-isoquant 0")
    val(modmasses)   // Comma-separated modification masses for --modlist (empty string = skip)
    path ionquant_dir // Unzipped IonQuant tool directory (optional). Pass [] when not using.

    output:
    tuple val(meta), path("${prefix}")                    , emit: results_dir
    tuple val(meta), path("${prefix}/ion.tsv")            , emit: ions, optional: true
    tuple val(meta), path("${prefix}/combined_ion.tsv")   , emit: combined_ions, optional: true
    tuple val(meta), path("${prefix}/combined_protein.tsv"), emit: combined_protein, optional: true
    tuple val(meta), path("${prefix}/combined_peptide.tsv"), emit: combined_peptide, optional: true
    tuple val(meta), path("${prefix}/combined_modified_peptide.tsv"), emit: combined_modified_peptide, optional: true
    tuple val(meta), path("${prefix}/combined_site_*.tsv") , emit: combined_site, optional: true
    tuple val(meta), path("${prefix}/*.tsv")              , emit: all_tsv, optional: true
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('ionquant'), eval("cat .ionquant_version"), emit: versions_ionquant, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def modlist_arg = modmasses ? '--modlist modmasses_ionquant.txt' : ''
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using this tool." >&2
        echo "Set ext.agree_fragpipe_license_agreement = true in your Nextflow config." >&2
        exit 1
    fi
    echo "INFO: User has set agree_fragpipe_license_agreement = true. Proceeding with licensed tool." >&2
    echo "User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)" > I_AGREE_FRAGPIPE_LICENSE_AGREEMENT

    TOOLS_DIR="${tools_dir}"

    # JAR discovery: prefer provided dir, fall back to tools_dir
    if [ -d "${ionquant_dir}" ] && [ "\$(ls -A ${ionquant_dir} 2>/dev/null)" ]; then
        IONQUANT_JAR=\$(find ${ionquant_dir} -maxdepth 1 -name "IonQuant*.jar" -type f 2>/dev/null | head -1 || true)
    else
        IONQUANT_JAR=\$(find "\$TOOLS_DIR" -name 'IonQuant*.jar' -type f 2>/dev/null | head -1 || true)
    fi
    JFREECHART_JAR=\$(find ${tools_dir} -name "jfreechart*.jar" -type f | sort | tail -1)

    # Version capture
    _iq_ver=\$(java -jar "\$IONQUANT_JAR" --version 2>&1 | grep -oP 'IonQuant-\\K[\\d.]+' | head -1 || true)

    echo "\$_iq_ver" > .ionquant_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version
    CP="\$IONQUANT_JAR"
    [ -n "\$JFREECHART_JAR" ] && CP="\$CP:\$JFREECHART_JAR"
    BRUKER_DIR="${tools_dir}/../ext/bruker"
    THERMO_DIR="${tools_dir}/../ext/thermo"
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.bruker.dir=\$BRUKER_DIR"
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.thermo.dir=\$THERMO_DIR"

    # Guard: when both quantification modes are disabled, copy dirs and exit cleanly.
    # This happens for TMT isobaric pass when no annotation file is provided —
    # IonQuant crashes with NullPointerException if run with both modes off.
    GUARD_CLI="${config_cli}"
    if echo "\$GUARD_CLI" | grep -q -- '--perform-isoquant 0' && echo "\$GUARD_CLI" | grep -q -- '--perform-ms1quant 0'; then
        echo "INFO: Both --perform-isoquant 0 and --perform-ms1quant 0 set. Skipping IonQuant, copying inputs to output." >&2
        for dir in */; do
            if [[ -f "\${dir}psm.tsv" && "\$dir" != "${prefix}/" && "\$dir" != "spectra/" ]]; then
                cp -rL "\$dir" "${prefix}/"
            fi
        done
        exit 0
    fi

    # Generate modlist file for MBR modification mass tracking (FragpipeRun.java:1852-1855)
    if [ -n "${modmasses}" ]; then
        echo "${modmasses}" | tr ',' '\\n' > modmasses_ionquant.txt
    fi

    # Generate filelist for IonQuant
    # Directories are already named by sample (from PHILOSOPHER_FILTER output)
    cat > filelist_ionquant.txt << EOF
flag	value
EOF

    # Find all psm.tsv files in sample directories
    for psm_file in */psm.tsv; do
        if [[ -f "\$psm_file" ]]; then
            echo -e "--psm\\t\$(pwd)/\$psm_file" >> filelist_ionquant.txt
        fi
    done

    echo -e "--specdir\\t\$(pwd)/spectra" >> filelist_ionquant.txt

    # Check if isoquant is enabled in config and annotation file is provided.
    # If isoquant is requested but no annotation file exists, disable it to avoid IonQuant failure.
    ISOQUANT_CLI="${config_cli}"
    if echo "\$ISOQUANT_CLI" | grep -q -- '--perform-isoquant 1'; then
        if [[ -n "${annotation_file}" && -f "${annotation_file}" ]]; then
            FIRST_LINE=\$(head -1 "${annotation_file}")
            if echo "\$FIRST_LINE" | grep -qP '^plex\\t'; then
                # Multi-plex TMTIntegrator format (6+ columns): extract per-plex 2-column annotations
                # Format: plex\tchannel\tsample\tsample_name\tcondition\treplicate
                # IonQuant needs 2-column: channel sample_name (per plex)
                for psm_file in */psm.tsv; do
                    if [[ -f "\$psm_file" ]]; then
                        PLEX_DIR=\$(dirname "\$psm_file")
                        ANNOT_FILE="\${PLEX_DIR}_annotation.txt"
                        # Extract channel and sample columns for this plex
                        awk -F'\\t' -v plex="\$PLEX_DIR" 'NR>1 && \$1==plex {print \$2 " " \$3}' "${annotation_file}" > "\$ANNOT_FILE"
                        if [[ -s "\$ANNOT_FILE" ]]; then
                            echo -e "--annotation\\t\$(pwd)/\$psm_file=\$(pwd)/\$ANNOT_FILE" >> filelist_ionquant.txt
                        fi
                    fi
                done
            else
                # Single-plex 2-column format: same annotation for all psm.tsv
                cp "${annotation_file}" annotation.txt
                for psm_file in */psm.tsv; do
                    if [[ -f "\$psm_file" ]]; then
                        echo -e "--annotation\\t\$(pwd)/\$psm_file=\$(pwd)/annotation.txt" >> filelist_ionquant.txt
                    fi
                done
            fi
        else
            # Isoquant requested but no annotation file provided — disable to avoid failure
            echo "WARNING: --perform-isoquant 1 requested but no annotation file provided. Disabling isoquant." >&2
            ISOQUANT_CLI=\$(echo "\$ISOQUANT_CLI" | sed 's/--perform-isoquant 1/--perform-isoquant 0/')
        fi
    fi

    export JAVA_OPTS="-Xmx${mem}G"
    CMD="java \$JAVA_OPTS \$NATIVE_FLAGS -cp \\"\$CP\\" ionquant.IonQuant --threads ${task.cpus} \$ISOQUANT_CLI ${modlist_arg} --multidir . --filelist filelist_ionquant.txt ${args}"
    rm -f ${prefix}/ionquant.log
    printf '%s\\n' "\$CMD" > ${prefix}/ionquant.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/ionquant.log

    # Copy combined TSVs (combined_ion.tsv, combined_protein.tsv, etc.) to output
    # Use cp instead of mv for Fusion filesystem compatibility
    shopt -s nullglob
    for f in *.tsv; do cp -L "\$f" ${prefix}/; done
    shopt -u nullglob

    # Copy IonQuant-updated per-sample directories to output.
    # IonQuant updates psm.tsv in-place, adding Intensity and tracing columns
    # (Apex Retention Time, Traced Scans, etc.) — matching FragPipe behavior
    # where IonQuant runs after philosopher report and modifies per-sample TSVs.
    # Use cp -rL instead of mv for Fusion filesystem compatibility (mv fails on Fusion symlinks).
    for dir in */; do
        if [[ -f "\${dir}psm.tsv" && "\$dir" != "${prefix}/" && "\$dir" != "spectra/" ]]; then
            cp -rL "\$dir" "${prefix}/"
        fi
    done

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    # Guard: when both quantification modes are disabled, copy dirs and exit cleanly.
    GUARD_CLI="${config_cli}"
    if echo "\$GUARD_CLI" | grep -q -- '--perform-isoquant 0' && echo "\$GUARD_CLI" | grep -q -- '--perform-ms1quant 0'; then
        for dir in */; do
            if [[ -f "\${dir}psm.tsv" && "\$dir" != "${prefix}/" && "\$dir" != "spectra/" ]]; then
                cp -rL "\$dir" "${prefix}/"
            fi
        done
        exit 0
    fi

    # Stub: copy staged per-sample directories to output (IonQuant would update these in-place)
    for dir in */; do
        if [[ -f "\${dir}psm.tsv" && "\$dir" != "${prefix}/" && "\$dir" != "spectra/" ]]; then
            cp -rL "\$dir" "${prefix}/"
        fi
    done

    touch ${prefix}/ion.tsv
    touch ${prefix}/combined_ion.tsv
    touch ${prefix}/combined_protein.tsv
    touch ${prefix}/combined_peptide.tsv
    touch ${prefix}/combined_modified_peptide.tsv
    touch ${prefix}/${prefix}_ionquant.log
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .ionquant_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}

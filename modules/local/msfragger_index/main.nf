/*
 * MSFRAGGER_INDEX: Prebuild MSFragger pepindex from FASTA + params (no spectra)
 *
 * Replicates CmdMsfraggerDigest.java from FragPipe: runs MSFragger with a params
 * file and FASTA but NO spectra files, which triggers digest-only mode and creates
 * pepindex files alongside the FASTA. These pepindex files can be reused by
 * downstream MSFRAGGER search tasks to skip redundant digest computation.
 *
 * Input:
 *   - meta: Sample/experiment metadata map
 *   - fasta: Protein sequence database (FASTA)
 *   - params_file: MSFragger params file (native fragger.params format)
 *
 * Output:
 *   - indexed_fasta: FASTA + co-located pepindex files (must stay together for MSFragger reuse)
 *   - versions: Software versions
 */
process MSFRAGGER_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)
    path params_file
    path msfragger_dir // Unzipped MSFragger tool directory (optional). Pass [] when not using.

    output:
    tuple val(meta), path(fasta), path("*.pepindex"), emit: indexed_fasta
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // CRITICAL: No 'def' for prefix — must be visible in output block (Nextflow 25.x scoping)
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def mem = task.ext.java_xmx ? task.ext.java_xmx as int : (task.memory ? task.memory.toGiga() : 8)
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false

    """
    export HOME=\$(pwd)
    export JAVA_OPTS="-Xmx${mem}G"

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using MSFragger." >&2
        echo "Set ext.agree_fragpipe_license_agreement = true in your Nextflow config." >&2
        exit 1
    fi
    echo "INFO: User has set agree_fragpipe_license_agreement = true. Proceeding with licensed tool." >&2
    echo "User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)" > I_AGREE_FRAGPIPE_LICENSE_AGREEMENT

    TOOLS_DIR="${tools_dir}"

    # JAR + native lib discovery: prefer provided dir, fall back to tools_dir
    if [ -d "${msfragger_dir}" ] && [ "\$(ls -A ${msfragger_dir} 2>/dev/null)" ]; then
        MSFRAGGER_JAR=\$(find ${msfragger_dir} -maxdepth 1 -name "MSFragger*.jar" -not -name "*-sources*" -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${msfragger_dir}/ext/bruker"
        THERMO_DIR="${msfragger_dir}/ext/thermo"
    else
        MSFRAGGER_JAR=\$(find "\$TOOLS_DIR" -name 'MSFragger*.jar' -not -name '*-sources*' -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${tools_dir}/../ext/bruker"
        THERMO_DIR="${tools_dir}/../ext/thermo"
    fi
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.bruker.dir=\$BRUKER_DIR"
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.thermo.dir=\$THERMO_DIR"

    # Version capture
    _ms_ver=\$(java -jar "\$MSFRAGGER_JAR" --help 2>&1 | grep -oP 'MSFragger-\\K[\\d.]+' | head -1 || true)

    echo "\$_ms_ver" > .msfragger_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Create local copy of params with correct database_name and calibrate_mass=0
    if [ -s "${params_file}" ]; then
        cp "${params_file}" fragger_index.params
        sed -i '/^database_name/d' fragger_index.params
        sed -i '/^num_threads/d' fragger_index.params
        sed -i '/^calibrate_mass/d' fragger_index.params
        # Ensure file ends with newline before appending
        sed -i -e '\$a\\' fragger_index.params
    else
        # No params file provided — create minimal with runtime params only
        touch fragger_index.params
    fi
    echo "database_name = \$(pwd)/${fasta}" >> fragger_index.params
    echo "num_threads = ${task.cpus}" >> fragger_index.params
    echo "calibrate_mass = 0" >> fragger_index.params

    # Run MSFragger with params + FASTA only (NO spectra files) -> digest-only mode
    # This creates *.pepindex file(s) alongside the FASTA
    java \$JAVA_OPTS -Dfile.encoding=UTF-8 \$NATIVE_FLAGS -jar "\$MSFRAGGER_JAR" fragger_index.params ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch "${fasta}.pepindex"
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .msfragger_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}

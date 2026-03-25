/*
 * MERGE_SPLIT_SEARCH: Merge MSFragger split-database search results per sample.
 *
 * For each sample, collects pepXML, PIN, and score histogram files from all
 * database chunks, sums histograms, generates expect functions via MSFragger,
 * and merges search results by re-ranking hits per spectrum.
 *
 * Faithfully ports the merge logic from FragPipe's msfragger_pep_split.py.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - chunk_pepxmls: pepXML files from all chunks for this sample
 *   - chunk_pins: PIN files from all chunks for this sample
 *   - chunk_histograms: Score histogram TSV files from all chunks
 *   - num_chunks: Number of database chunks
 *   - params_file: MSFragger params file (for output_report_topN, output_max_expect)
 *
 * Output:
 *   - pepxml: Merged pepXML file
 *   - pin: Merged PIN file
 *   - versions: Software versions
 */
process MERGE_SPLIT_SEARCH {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(chunk_pepxmls), path(chunk_pins), path(chunk_histograms)
    val num_chunks
    path params_file
    path msfragger_dir // Unzipped MSFragger tool directory (optional). Pass [] when not using.

    output:
    tuple val(meta), path("merged/*.pepXML"), emit: pepxml
    tuple val(meta), path("merged/*.pin")   , emit: pin
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_report_topN = task.ext.output_report_topN ?: 1
    def output_max_expect = task.ext.output_max_expect ?: 50.0
    def fasta_path = task.ext.fasta_path ?: ''
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false
    def mem = task.ext.java_xmx ? task.ext.java_xmx as int : (task.memory ? task.memory.toGiga() : 8)

    // Organize chunk-prefixed files into subdirectories for merge_split_search.py.
    // MSFRAGGER outputs are prefixed as chunk_{N}_{sample}.pepXML etc. to avoid
    // Nextflow staging collisions. Here we strip the prefix and sort into chunk_N/ dirs.
    """
    export JAVA_OPTS="-Xmx${mem}G"

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using this tool." >&2
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

    MSFRAGGER_CMD="java \$JAVA_OPTS -Dfile.encoding=UTF-8 \$NATIVE_FLAGS -jar \$MSFRAGGER_JAR"

    # Create chunk directories
    for i in \$(seq 0 \$((${num_chunks} - 1))); do
        mkdir -p "chunk_\${i}"
    done

    # Sort chunk-prefixed files into chunk directories, stripping the prefix.
    # Files are named: chunk_{N}_{original_name}
    for f in ${chunk_pepxmls} ${chunk_pins} ${chunk_histograms}; do
        basename_f=\$(basename "\$f")
        if [[ "\$basename_f" =~ ^chunk_([0-9]+)_(.*)\$ ]]; then
            chunk_idx="\${BASH_REMATCH[1]}"
            original_name="\${BASH_REMATCH[2]}"
            cp "\$f" "chunk_\${chunk_idx}/\${original_name}"
        else
            echo "WARNING: File '\$basename_f' does not have chunk prefix, skipping" >&2
        fi
    done

    # Build comma-separated chunk directory list
    CHUNK_DIRS=\$(printf "chunk_%d," \$(seq 0 \$((${num_chunks} - 1))) | sed 's/,\$//')

    # Run merge
    merge_split_search.py \\
        --sample_name ${prefix} \\
        --chunk_dirs \${CHUNK_DIRS} \\
        --num_chunks ${num_chunks} \\
        --msfragger_cmd "\$MSFRAGGER_CMD" \\
        --outdir merged \\
        --output_report_topN ${output_report_topN} \\
        --output_max_expect ${output_max_expect} \\
        --fasta_path '${fasta_path}'

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p merged
    touch "merged/${prefix}.pepXML"
    touch "merged/${prefix}.pin"
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .msfragger_version
    echo "\$_fp_ver" > .fragpipe_version
    """
}

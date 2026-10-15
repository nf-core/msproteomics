/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FragPipe module input helpers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pure functions that build the exact paths, argv words and file bodies the FragPipe
    modules run with, so that a module script only runs its tool.
    Modules include them from their own main.nf:
        include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/**
 * Paths of the FragPipe 24.0 tools under a FragPipe tools directory.
 *
 * The names are the ones the FragPipe 24.0 release runs (the argv of the golden FragPipe runs), under
 * task.ext.fragpipe_tools_dir (nextflow.config); the MSFragger, IonQuant and diaTracer jars, whose names
 * differ between distributions, are found in the task (fragpipeJarTokens, fragpipeJarLines). FragPipe passes the
 * license it finds in <tools>/.. (FragpipeLocations.locateLicense); MSFragger, IonQuant, MBG,
 * PTM-Shepherd and diaTracer get the Thermo/Bruker libraries from <MSFragger dir>/ext
 * (FragpipeRun.java:1348-1356). Nothing here touches the file system: the tools run in the container,
 * where a missing file fails the tool.
 *
 * @param task Task whose ext.fragpipe_tools_dir (set in nextflow.config) is the FragPipe tools directory,
 *             e.g. /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools
 * @return Map of tool name to absolute path (String), plus fragpipe_version and the tool versions
 */
def fragpipeTools(task) {
    def tools_dir = task.ext.fragpipe_tools_dir
    def tools = (tools_dir ?: '').toString().replaceAll(/\/+$/, '')
    if (!tools.startsWith('/')) {
        error("fragpipeTools: task.ext.fragpipe_tools_dir must be the absolute FragPipe tools directory, got '${tools_dir}'")
    }
    def root = tools.substring(0, tools.lastIndexOf('/'))
    def version_match = tools =~ /fragpipe-(\d+\.\d+[0-9.]*)/
    if (!version_match.find()) {
        error("fragpipeTools: task.ext.fragpipe_tools_dir '${tools_dir}' names no fragpipe-<version> directory; the FragPipe jar is lib/fragpipe-<version>.jar (e.g. /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools)")
    }
    def fragpipe_version = version_match.group(1)
    def versions = [crystalc: '1.5.10', diaumpire: '2.3.3', mbg: '0.3.6', fpmeta: '1.0.1', fragannotator: '1.0.0']
    return [
        tools_dir        : tools,
        // the tools path FragPipe writes into argv and config files: getDirTools() = <install>/lib/../tools,
        // not normalised (FragpipeLocations.java:105; CmdMBGMatch.java:159, CmdOPair.java:148, PTMSGlycanAssignPanel.java:131)
        argv_tools_dir   : "${root}/lib/../tools".toString(),
        fragpipe_root    : root,
        fragpipe_version : fragpipe_version,
        fragpipe_lib     : "${root}/lib".toString(),
        fragpipe_jar     : "${root}/lib/fragpipe-${fragpipe_version}.jar".toString(),
        license          : "${root}/license.dat".toString(),
        // FragpipeLocations.getDirDatabases(): <install>/lib/../databases, normalised (FragpipeLocations.java:65, :109)
        databases_dir    : "${root}/databases".toString(),
        // MSFragger, IonQuant and diaTracer jars are named by the distribution (MSFragger-4.4.1.jar, or with an
        // edition word before the version): the task finds them as FragPipe does (fragpipeJarLines), so these
        // are tokens that shellArgs, heredocText, shellQuote and taskFileCommand turn into its shell variables
        // *_version: the FragPipe 24.0 release versions of those jars, for stubs (a run writes the found jar's version)
        msfragger_jar    : '{MSFRAGGER_JAR}',
        msfragger_version: '4.4.1',
        thermo_dir       : '{MSFRAGGER_DIR}/ext/thermo',
        bruker_dir       : '{MSFRAGGER_DIR}/ext/bruker',
        ionquant_jar     : '{IONQUANT_JAR}',
        ionquant_version : '1.11.20',
        diatracer_jar    : '{DIATRACER_JAR}',
        diatracer_version: '2.2.1',
        jfreechart_jar   : "${tools}/jfreechart-1.5.3.jar".toString(),
        // batmass-io and MSBooster: the versions the image's FragPipe build names (ToolingUtils.BATMASS_IO_JAR,
        // CmdMSBooster.JAR_MSBOOSTER_NAME), read in the task (fragpipeJarLines)
        batmass_io_jar   : '{BATMASS_IO_JAR}',
        grppr_jar        : "${tools}/grppr-0.3.23.jar".toString(),
        // Versions FragPipe 24.0 compiles in and names the jars by (it reads no version from a jar):
        // CmdCrystalc.CRYSTALC_VERSION (CmdCrystalc.java:48-49), UmpireParams.UMPIRESE_VERSION
        // (UmpireParams.java:80-81), CmdMBGMatch.JAR_MBG_NAME (CmdMBGMatch.java:50), CmdMetaproteomics.VERSION
        // (CmdMetaproteomics.java:33-34), TabRun.FRAG_ANNOTATOR (TabRun.java:96); the jars' pom.properties
        // (crystalc, FragPipe-Meta, FragAnnotator) say the same
        crystalc_jar     : "${tools}/original-crystalc-${versions.crystalc}.jar".toString(),
        crystalc_version : versions.crystalc,
        diaumpire_jar    : "${tools}/DIA_Umpire_SE-${versions.diaumpire}.jar".toString(),
        diaumpire_version: versions.diaumpire,
        mbg_jar          : "${tools}/MBG-${versions.mbg}.jar".toString(),
        mbg_version      : versions.mbg,
        msbooster_jar    : '{MSBOOSTER_JAR}',
        msbooster_version: '1.4.17',
        diann            : "${tools}/diann/1.8.2_beta_8/linux/diann-1.8.1.8".toString(),
        diann_so_dir     : "${tools}/diann_so".toString(),
        unimod_obo       : "${tools}/unimod.obo".toString(),
        glycan_residues  : "${tools}/Glycan_Databases/glycan_residues.txt".toString(),
        glycan_mods      : "${tools}/Glycan_Databases/glycan_mods.txt".toString(),
        fragreporter_jar : "${tools}/FragReporter-1.1.2.jar".toString(),
        fragannotator_jar    : "${tools}/FP-PDV/FragAnnotator-${versions.fragannotator}.jar".toString(),
        fragannotator_version: versions.fragannotator,
        fpmeta_jar       : "${tools}/metaproteomics/FP-Meta-${versions.fpmeta}.jar".toString(),
        fpmeta_version   : versions.fpmeta,
        taxon_names      : "${tools}/metaproteomics/names.dmp".toString(),
        taxon_nodes      : "${tools}/metaproteomics/nodes.dmp".toString(),
        fpop_script      : "${tools}/fpop/FragPipe_FPOP_Analysis.py".toString(),
        philosopher      : "${tools}/Philosopher/philosopher-v5.1.3-RC9".toString(),
    ]
}

/**
 * The fragpipeTools tokens of the distribution-named jars and the task shell variable each stands for.
 *
 * @return Map of token to shell variable name
 */
def fragpipeJarTokens() {
    return ['{MSFRAGGER_JAR}': 'MSFRAGGER_JAR', '{MSFRAGGER_DIR}': 'MSFRAGGER_DIR', '{IONQUANT_JAR}': 'IONQUANT_JAR', '{DIATRACER_JAR}': 'DIATRACER_JAR',
            '{BATMASS_IO_JAR}': 'BATMASS_IO_JAR', '{MSBOOSTER_JAR}': 'MSBOOSTER_JAR', '{FRAGPIPE_HASHSET_CP}': 'FRAGPIPE_HASHSET_CP']
}

/**
 * Text with every fragpipeJarTokens token as its shell variable (${MSFRAGGER_JAR}, ...), for text that bash
 * expands (double-quoted words, unquoted here-documents); in a single-quoted word, open and close end the
 * quoting around the variable ("'\"" and "\"'").
 *
 * @param text  text already escaped for bash
 * @param open  text before each variable
 * @param close text after each variable
 * @return text with the tokens replaced
 */
def expandJarTokens(text, open = '', close = '') {
    def out = text.toString()
    fragpipeJarTokens().each { token, variable -> out = out.replace(token, "${open}\${${variable}}${close}") }
    return out
}

/**
 * Shell lines that find the MSFragger, IonQuant and diaTracer jars in the FragPipe tools folder as FragPipe
 * does (TabConfig.getJarPath, TabConfig.java:176-178, :188-209): every file under the folder whose name
 * matches <tool>-(.*).jar, case-insensitive, the highest version winning; MSFragger's ext/thermo and
 * ext/bruker are next to its jar (FragpipeRun.java:1348-1356). So MSFragger-4.4.1.jar (academic) and an
 * edition-named jar (FragPipePlus) both work. Sets MSFRAGGER_JAR, MSFRAGGER_DIR, IONQUANT_JAR and
 * DIATRACER_JAR and writes .<tool>_version from the jar name. batmass-io and MSBooster are the file names the
 * image's FragPipe build has compiled in (ToolingUtils.BATMASS_IO_JAR, CmdMSBooster.JAR_MSBOOSTER_NAME; the
 * FragPipe 24.0 release names batmass-io-1.36.5 / MSBooster-1.4.14, FragPipePlus 1.36.6 / 1.4.17), read from
 * its lib/fragpipe-<version>.jar, under the tools folder (FragpipeLocations.checkToolsMissing): BATMASS_IO_JAR,
 * MSBOOSTER_JAR. A required tool without a jar stops the task with an error naming the folder; an optional
 * one leaves its variables empty.
 *
 * Approved exception to "modules only run tools" (docs/CONTRIBUTING.md, Pipeline specific contribution
 * guidelines): the licensed FragPipe image is a runtime-only parameter (--fragpipe_container), so these jar
 * names can only be resolved inside the task. The lookup copies FragPipe's (TabConfig.java:176-178, :188-209;
 * ToolingUtils.java:62; CmdMSBooster.java:52; FragpipeRun.java:1348-1356); it is not a pattern for other
 * in-task lookups.
 *
 * @param tools    fragpipeTools() result
 * @param required tools the task runs: any of 'msfragger', 'ionquant', 'diatracer', 'batmass_io', 'msbooster'
 * @param optional tools the task uses only when present (e.g. MSFragger's Thermo libraries)
 * @return shell lines
 */
def fragpipeJarLines(tools, required, optional = []) {
    def names = [msfragger: 'MSFragger', ionquant: 'IonQuant', diatracer: 'diaTracer', batmass_io: 'batmass-io', msbooster: 'MSBooster']
    // the FragPipe class holding each build-named jar, and the jar name pattern in its constant pool
    def constants = [batmass_io: ['org/nesvilab/fragpipe/cmd/ToolingUtils.class', 'batmass-io-[0-9.]+[.]jar'],
                     msbooster : ['org/nesvilab/fragpipe/cmd/CmdMSBooster.class', 'MSBooster-[0-9.]+[.]jar']]
    def unknown = (required + optional) - names.keySet()
    if (unknown) {
        error("fragpipeJarLines: unknown tool(s) ${unknown}; known: ${names.keySet().join(', ')}")
    }
    def dir = shellArgs([tools.tools_dir])
    def lines = [
        '# FragPipe\'s jar lookup in its tools folder (TabConfig.java:188-209): the highest <tool>-<version>.jar',
        // a missing tools folder finds nothing (the required-jar check names it): find's status must not stop a pipefail shell
        "fragpipe_jar() { { find ${dir} -type f -iname \"\$1-*.jar\" 2>/dev/null || true; } | awk -F/ '{ print \$NF \"\\t\" \$0 }' | sort -V -k1,1 | tail -n 1 | cut -f 2; }",
    ]
    if ((required + optional).any { key -> key in constants }) {
        lines << "fragpipe_build_jar() { local n; n=\$({ unzip -p ${shellArgs([tools.fragpipe_jar])} \"\$1\" 2>/dev/null || true; } | { grep -aoE \"\$2\" || true; } | sort -u | tail -n 1); [ -n \"\$n\" ] && [ -f ${dir}/\"\$n\" ] && echo ${dir}/\"\$n\" || true; }"
    }
    (required + optional).unique().each { key ->
        def variable = "${key.toUpperCase()}_JAR"
        lines << (key in constants ? "${variable}=\$(fragpipe_build_jar ${constants[key][0]} '${constants[key][1]}')" : "${variable}=\$(fragpipe_jar ${key})")
        if (key in required) {
            lines << "[ -n \"\$${variable}\" ] || { echo \"ERROR: no ${names[key]}-*.jar under ${tools.tools_dir} (ext.fragpipe_tools_dir; run in the FragPipe image given with --fragpipe_container)\" >&2; exit 1; }"
            if (key != 'batmass_io') {
                lines << "basename \"\$${variable}\" .jar | grep -oE '[0-9]+([.][0-9]+)+\$' > .${key}_version"
            }
        }
        if (key == 'msfragger') {
            lines << 'MSFRAGGER_DIR=${MSFRAGGER_JAR:+$(dirname "$MSFRAGGER_JAR")}'
        }
    }
    return lines.join('\n')
}

/**
 * Shell lines that set FRAGPIPE_HASHSET_CP ({FRAGPIPE_HASHSET_CP}) to the class path FragPipe joins from a
 * java.util.HashSet of its words (default capacity 16): ascending bucket (h ^ (h >>> 16)) & 15 of
 * String.hashCode, insertion order within a bucket. The words are the task's paths (jar tokens expanded), so
 * the order is the one FragPipe gets for the same files (CmdWriteSubMzml.java:93-95).
 *
 * @param words class path entries in insertion order, may hold jar tokens
 * @return shell lines
 */
def hashSetClassPathLines(words) {
    return [
        "fp_bucket() { local s=\$1 h=0 i c; for ((i = 0; i < \${#s}; i++)); do printf -v c '%d' \"'\${s:i:1}\"; h=\$(( (31 * h + c) & 0xffffffff )); done; echo \$(( (h ^ (h >> 16)) & 15 )); }",
        "FP_CP_WORDS=(${shellArgs(words)})",
        'FRAGPIPE_HASHSET_CP=$(for i in "${!FP_CP_WORDS[@]}"; do printf \'%s\\t%s\\t%s\\n\' "$(fp_bucket "${FP_CP_WORDS[$i]}")" "$i" "${FP_CP_WORDS[$i]}"; done | sort -n -k1,1 -k2,2 | cut -f 3 | paste -sd :)',
    ].join('\n')
}

/**
 * Set values of "key = value" lines in a FragPipe params text, keeping each line's tab comment.
 *
 * FragPipe writes params as "key = value\t\t\t# comment" (Props.java:240-273) and fills database_name,
 * num_threads (CmdMsfragger.java:441, TabMsfragger.java:1657) or check_spectral_files
 * (msfragger_pep_split.py) in place. Every line that starts with "<key> = " gets the new value; the text
 * from the line's first tab on is kept. Lines are separated by '\n' and every line ends with '\n'.
 * A key without a line is an error: the params text is not the one FragPipe fills.
 *
 * @param text   params text
 * @param values Map of key to new value
 * @return the params text with the values replaced
 */
def setFragpipeParamValues(text, Map values) {
    def lines = text.toString().split('\n', -1) as List
    if (lines && lines.last() == '') {
        lines = lines.dropRight(1)
    }
    def missing = values.keySet().findAll { k -> !lines.any { line -> line.startsWith("${k} = ") } }
    if (missing) {
        error("setFragpipeParamValues: no '<key> = ' line for ${missing} in the params text")
    }
    def out = lines.collect { line ->
        def key = values.keySet().find { k -> line.startsWith("${k} = ") }
        if (key == null) {
            return line
        }
        def tab = line.indexOf('\t')
        return "${key} = ${values[key]}" + (tab >= 0 ? line.substring(tab) : '')
    }
    return out.collect { line -> "${line}\n" }.join('')
}

/**
 * Body of an unquoted bash here-document that writes `text` byte for byte, with the token {TASK_DIR}
 * standing for the task's work directory ($PWD in the task).
 *
 * The characters bash expands in an unquoted here-document (\, $ and `) are escaped, so the only
 * expansion left is ${PWD} ({DIANN_INPUT_DIR} becomes ${DIANN_INPUT_DIR}, see diannInputDirLines). One trailing newline is dropped because the here-document adds it back.
 *
 * @param text file content, may contain {TASK_DIR}
 * @return here-document body
 */
def heredocText(text) {
    def body = text.toString()
    if (body.endsWith('\n')) {
        body = body.substring(0, body.length() - 1)
    }
    return body
        .replace('\\', '\\\\')
        .replace('$', '\\$')
        .replace('`', '\\`')
        .replace('{TASK_DIR}', '${PWD}')
        .replace('{DIANN_INPUT_DIR}', '${DIANN_INPUT_DIR}')
        .with { t -> expandJarTokens(t) }
}

/**
 * Runtime keys FragPipe appends to the Crystal-C algorithm settings (CmdCrystalc.java:190-196).
 *
 * The text starts with an empty line, as the block the module appended before; paths are under the
 * task directory ({TASK_DIR}, see heredocText).
 *
 * @param threads  thread count (the task's CPUs)
 * @param fasta    staged protein database name
 * @param lcms     staged LC-MS file name; raw_file_extension is its extension, case kept
 *                 (StringUtils.afterLastDot, CmdCrystalc.java:192-193)
 * @param out_dir  output directory under the task directory
 * @return text appended to the Crystal-C params file
 */
def crystalcRuntimeParams(threads, fasta, lcms, out_dir) {
    def raw_ext = lcms.toString().tokenize('.').last()
    return [
        '',
        "thread = ${threads}",
        "fasta = {TASK_DIR}/${fasta}",
        'raw_file_location = {TASK_DIR}',
        "raw_file_extension = ${raw_ext}",
        "output_location = {TASK_DIR}/${out_dir}",
    ].collect { line -> "${line}\n" }.join('')
}

/**
 * Shell lines that replace staged links by copies of what they name (cp -rL, then the copy takes the link's
 * place), for a tool that misreads a link or writes through it. Only these inputs are copied; the others stay
 * links. A staged directory (a Waters .raw folder) is copied whole; rm + mv replaces the link itself, where
 * mv onto a link to a directory would move the copy into that directory. Names are quoted with shellArgs.
 *
 * @param names staged file or directory names, relative to the task directory (String, Path or a List of them)
 * @return one line per name, joined by newlines
 */
def linkCopies(names) {
    def list = names instanceof List ? names : [names]
    return list.collect { n -> n instanceof java.nio.file.Path ? n.name : n.toString() }.collect { n ->
        def f = shellArgs([n])
        def c = shellArgs(["${n}.copy"])
        "[ ! -L ${f} ] || { cp -rL ${f} ${c} && rm ${f} && mv ${c} ${f}; }"
    }.join('\n')
}

/**
 * linkCopies() of the staged inputs whose absolute path in `dir` contains `needle`, decided when the task
 * runs (the task directory is known only then). For Philosopher's PeptideProphet, which rewrites in place
 * every input whose filepath.Abs path contains its output name, "interact" with FragPipe's options
 * (philosopher lib/ext/peptideprophet/peptideprophet.go:80-96, :457): the other inputs stay links.
 *
 * @param names  staged file names, relative to the task directory (String, Path or a List of them)
 * @param dir    directory, relative to the task directory, the tool reads the inputs from ('.' = task dir)
 * @param needle text that makes the tool write the input in place
 * @return one line per name, joined by newlines
 */
def linkCopiesWhenPathContains(names, dir, needle) {
    def list = names instanceof List ? names : [names]
    return list.collect { n -> n instanceof java.nio.file.Path ? n.name : n.toString() }.collect { n ->
        "case \"\$PWD\"/${shellArgs(["${dir}/${n}"])} in *${shellArgs([needle])}*) ${linkCopies(n)} ;; esac"
    }.join('\n')
}

/**
 * linkCopies() of the Thermo .raw files among `names`: batmass-io reads .raw under Mono, where
 * FileInfo.Length of a symlink is the link size, so a staged link reads as an empty file.
 *
 * @param names staged file names (String, Path or a List of them)
 * @return one line per .raw file ('' when there is none)
 */
def rawLinkCopies(names) {
    def list = names instanceof List ? names : [names]
    return linkCopies(list.findAll { n -> (n instanceof java.nio.file.Path ? n.name : n.toString()).toLowerCase().endsWith('.raw') })
}

/**
 * Shell lines that write `text` to a file byte for byte (a here-document of heredocText(text)).
 *
 * A text that does not end with a newline is written without one (the here-document's newline is
 * truncated), as FragPipe's FileUtils.write does for the manifest. {TASK_DIR} becomes the task directory.
 *
 * @param name   file name, relative to the task directory
 * @param text   file content
 * @param append true to append instead of overwrite
 * @return shell lines
 */
def writeTextFile(name, text, append = false) {
    def body = text.toString()
    def redirect = append ? '>>' : '>'
    if (body.isEmpty()) {
        return append ? '' : ": > '${name}'"
    }
    def lines = ["cat ${redirect} '${name}' <<FP_EOF", heredocText(body), 'FP_EOF']
    if (!body.endsWith('\n')) {
        lines << "truncate -s -1 '${name}'"
    }
    return lines.join('\n')
}

/**
 * Text of FragPipe's fragpipe-files.fp-manifest: one "<path>\t<experiment>\t<bioreplicate>\t<data type>"
 * line per LC-MS file, no newline after the last line (TabWorkflow.java:1241-1250, FileUtils.write).
 *
 * @param rows        [name, experiment, bioreplicate, data type] per file; null cells are written empty
 * @param spectra_dir directory the files are staged in, e.g. '{TASK_DIR}/spectra'
 * @return manifest text
 */
def fragpipeManifestText(rows, spectra_dir) {
    return rows.collect { row ->
        def cells = row.collect { cell -> cell == null ? '' : cell.toString() }
        (["${spectra_dir}/${cells[0]}"] + cells.drop(1)).join('\t')
    }.join('\n')
}

/**
 * Replace the run-time tokens of a parser argv (e.g. @THREADS@, @TOOLS_DIR@) in every word.
 *
 * @param argv   argv words from parse_fragpipe_workflow.py
 * @param tokens Map of token to value
 * @return argv words with the tokens replaced
 */
def fillArgvTokens(argv, Map tokens) {
    return argv.collect { word ->
        def w = word.toString()
        tokens.each { token, value -> w = w.replace(token.toString(), value.toString()) }
        w
    }
}

/**
 * argv words as bash words: each word double-quoted with \, $, ` and " escaped, {TASK_DIR} as ${PWD}.
 *
 * @param words argv words
 * @return words joined by spaces, ready to place in a bash array or command line
 */
def shellArgs(words) {
    return words.collect { word ->
        def w = word.toString()
            .replace('\\', '\\\\')
            .replace('$', '\\$')
            .replace('`', '\\`')
            .replace('"', '\\"')
            .replace('{TASK_DIR}', '${PWD}')
        "\"${expandJarTokens(w)}\""
    }.join(' ')
}

/**
 * Fail unless the user agreed to the FragPipe license (task.ext.agree_fragpipe_license_agreement), as the
 * licensed tools (MSFragger, IonQuant, diaTracer) require.
 *
 * @param task Task whose ext.agree_fragpipe_license_agreement must be true
 * @return the text of the I_AGREE_FRAGPIPE_LICENSE_AGREEMENT record file
 */
def requireFragpipeLicenseAgreement(task) {
    if (task.ext.agree_fragpipe_license_agreement != true) {
        error("You must agree to the FragPipe license agreement before using this tool. Set ext.agree_fragpipe_license_agreement = true in your Nextflow config.")
    }
    return 'User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)\n'
}

/**
 * Value of "key = value" in a FragPipe params text: the text after "=", up to a tab or '#', trimmed.
 *
 * @param text params text
 * @param key  parameter name
 * @return the value, or null when the params have no such line
 */
def fragpipeParamValue(text, key) {
    def name = key.toString()
    def line = text.toString().readLines().find { l -> l.startsWith(name) && l.substring(name.length()).replaceAll(/^ */, '').startsWith('=') }
    if (line == null) {
        return null
    }
    return line.substring(line.indexOf('=') + 1).replaceAll(/[\t#].*$/, '').trim()
}

/**
 * Merge settings msfragger_pep_split.py reads from the search params (FragPipe 24.0, lines 76-80):
 * output_file_extension (default pepXML), output_report_topN (required) and output_max_expect
 * (default 50.0), each the first word after "=".
 *
 * @param params_text MSFragger params text of the file
 * @return Map with output_file_extension, output_report_topN, output_max_expect (Strings)
 */
def msfraggerSplitMergeSettings(params_text) {
    def words = ['output_file_extension', 'output_report_topN', 'output_max_expect'].collectEntries { key ->
        def value = fragpipeParamValue(params_text, key)
        [key, value ? value.tokenize()[0] : null]
    }
    def top_n = words.output_report_topN
    if (top_n == null || !(top_n ==~ /\d+/)) {
        error("msfragger_pep_split.py needs output_report_topN in the MSFragger params, found '${top_n}'")
    }
    return [
        output_file_extension: words.output_file_extension ?: 'pepXML',
        output_report_topN   : top_n,
        output_max_expect    : words.output_max_expect ?: '50.0',
    ]
}

/**
 * Files an MSFragger search writes, for the stub (CmdMsfragger.java:103-126): <stem>.pepXML for DDA
 * params or <stem>_rank1..N.pepXML (N = output_report_topN) otherwise, <stem>.pin (output_format with pin) and
 * <stem>_scores_histogram.tsv in out_dir; next to the LC-MS file <stem>_uncalibrated.mzML for .raw/.d
 * and <stem>_calibrated.mzML with write_calibrated_mzml = 1 and calibration on, <stem>.mzBIN_all with
 * write_mzbin_all = 1; one
 * <fasta>.<N>.pepindex per digestion set (2 with calibrate_mass 1/2).
 *
 * @param spectra     LC-MS file names
 * @param fasta       FASTA file name
 * @param params_text MSFragger params text
 * @param out_dir     directory of the search outputs
 * @return relative paths of the files
 */
def msfraggerStubFiles(spectra, fasta, params_text, out_dir) {
    def values = ['data_type', 'output_report_topN', 'calibrate_mass', 'write_calibrated_mzml', 'output_format'].collectEntries { key ->
        def v = fragpipeParamValue(params_text, key)
        if (v == null) {
            error("MSFragger params have no ${key}")
        }
        [key, v]
    }
    def data_type = values.data_type
    def top_n = values.output_report_topN as int
    def calibrate_mass = values.calibrate_mass as int
    def write_calibrated = values.write_calibrated_mzml
    def names = (spectra instanceof List ? spectra : [spectra]).collect { f -> f instanceof java.nio.file.Path ? f.name : f.toString() }
    def files = names.collectMany { name ->
        def stem = name.contains('.') ? name.substring(0, name.lastIndexOf('.')) : name
        def pepxml = data_type != '0' ? (1..top_n).collect { r -> "${out_dir}/${stem}_rank${r}.pepXML".toString() } : ["${out_dir}/${stem}.pepXML".toString()]
        def lcms = []
        if (name.toLowerCase() ==~ /.*\.(d|raw)$/) {
            lcms << "${stem}_uncalibrated.mzML".toString()
        }
        if (write_calibrated == '1' && calibrate_mass > 0) {
            lcms << "${stem}_calibrated.mzML".toString()
        }
        // write_mzbin_all = 1: <stem>.mzBIN_all next to the LC-MS file (WriteSubMzml.java:73)
        if (fragpipeParamValue(params_text, 'write_mzbin_all') == '1') {
            lcms << "${stem}.mzBIN_all".toString()
        }
        // MSFragger writes <stem>.pin only when output_format names it (tsv_pepXML writes none)
        def pin = values.output_format.contains('pin') ? ["${out_dir}/${stem}.pin".toString()] : []
        pepxml + pin + ["${out_dir}/${stem}_scores_histogram.tsv".toString()] + lcms
    }
    return files + ["${fasta}.1.pepindex".toString()] + (calibrate_mass != 0 ? ["${fasta}.2.pepindex".toString()] : [])
}

/**
 * MSFragger options FragPipe appends after the LC-MS files (CmdMsfragger.java:586, :628-630):
 *
 *     boolean addNoPtmScoreInPin = tabRun.isWriteSubMzml();
 *     ...
 *     if (addNoPtmScoreInPin) { cmd.add("--no-ptm-score-in-pin"); }
 *
 * FragPipe sets write_mzbin_all in every params set from the same checkbox (FragpipeRun.java:1481,
 * CmdMsfragger.java:443 params.setWriteMzbinAll(writeMzbinAll)), so write_mzbin_all = 1 in the params
 * text means Write sub mzML is on. In the second pass FragPipe adds the flag when an LC-MS file of the
 * command is a _sub.mzML with its _percolator.weights below its directory (CmdMsfragger.java:612-626;
 * subMzmlWeights sets meta.no_ptm_score_in_pin).
 *
 * @param params_text         MSFragger params text of the file
 * @param meta                meta of the file (meta.no_ptm_score_in_pin: second pass)
 * @return [] or ['--no-ptm-score-in-pin']
 */
def msfraggerTrailingArgv(params_text, Map meta) {
    return (fragpipeParamValue(params_text, 'write_mzbin_all') == '1' || meta.no_ptm_score_in_pin == true) ? ['--no-ptm-score-in-pin'] : []
}

/**
 * WriteSubMzml command as CmdWriteSubMzml.configure builds it (CmdWriteSubMzml.java:82-130):
 *
 *     java -Xmx<ram>G [-Dlibs.bruker.dir=<dir>] [-Dlibs.thermo.dir=<dir>] -cp <classpath>
 *         org.nesvilab.fragpipe.util.WriteSubMzml <LC-MS file> <wd>/<group>/psm.tsv <wd>/<stem>_sub.mzML <threshold> 1
 *
 * run in <wd>/<group> (pb.directory). The Bruker and Thermo libraries are MSFragger's ext/bruker and
 * ext/thermo (FragpipeRun.java:1352-1353). The class path joins a HashSet of the batmass-io jar and
 * <install>/lib/* (:93-95): the task orders it as that java.util.HashSet does (hashSetClassPathLines, classpath_words). The threshold
 * is TabRun.getSubMzmlProbThreshold(), a float printed with probabilityThreshold + "" (:126; TabRun.java:495-496:
 * ((Double) spinner.getValue()).floatValue()). The trailing 1 deletes the .mzBIN_all on exit (WriteSubMzml.java:64, :75-77).
 * -Duser.home: the container user has no passwd entry.
 * Paths are under {TASK_DIR}, the task directory standing for FragPipe's work dir; the LC-MS file and its
 * .mzBIN_all are staged in lcms_dir.
 *
 * @param lcms_name      LC-MS file name, staged in lcms_dir
 * @param lcms_dir       directory of the LC-MS file under the task dir
 * @param meta           meta of the LC-MS file (meta.group: FragPipe group, '' for the blank group: psm.tsv in the work dir)
 * @param prob_threshold workflow text of tab-run.sub_mzml_prob_threshold (tool_configs writesubmzml.prob_threshold)
 * @param ram_gb         Java heap in GB
 * @param tools          fragpipeTools()
 * @return Map: argv (words), classpath_words, group_dir (cwd relative to the task dir), psm (relative path of psm.tsv), sub_mzml (output name)
 */
def writeSubMzmlArgv(lcms_name, lcms_dir, Map meta, prob_threshold, ram_gb, Map tools) {
    def group = meta.group
    if (group == null) {
        error("writeSubMzmlArgv: the FragPipe group of ${lcms_name} is missing")
    }
    def threshold = Float.toString(Double.parseDouble(prob_threshold.toString()).floatValue())
    // the HashSet order of [batmass-io jar, <install>/lib/*] depends on the image's batmass-io name: the task
    // computes it (hashSetClassPathLines(classpath_words))
    def classpath_words = [tools.batmass_io_jar.toString(), "${tools.fragpipe_lib}/*".toString()]
    def name = lcms_name.toString()
    def stem = name.lastIndexOf('.') < 0 ? name : name.substring(0, name.lastIndexOf('.'))
    def group_dir = group.toString() == '' ? '.' : group.toString()
    def psm = group_dir == '.' ? 'psm.tsv' : "${group_dir}/psm.tsv".toString()
    def sub_mzml = "${stem}_sub.mzML".toString()
    def argv = ['java', '-Duser.home={TASK_DIR}', "-Xmx${ram_gb}G",
        "-Dlibs.bruker.dir=${tools.bruker_dir}", "-Dlibs.thermo.dir=${tools.thermo_dir}",
        '-cp', '{FRAGPIPE_HASHSET_CP}', 'org.nesvilab.fragpipe.util.WriteSubMzml',
        "{TASK_DIR}/${lcms_dir}/${name}", "{TASK_DIR}/${psm}", "{TASK_DIR}/${sub_mzml}", threshold, '1']
    return [argv: argv.collect { w -> w.toString() }, classpath_words: classpath_words, group_dir: group_dir, psm: psm, sub_mzml: sub_mzml]
}

/**
 * MSBooster's thread count: the numThreads line parse_fragpipe_workflow.py writes into the MSBooster params
 * (workflow.threads / --threads; MSBooster's RT loess output depends on it, CmdMSBooster.java:164). It sizes the
 * MSBOOSTER task and is checked against the params again when the file is rendered (msboosterParamsText).
 *
 * @param params_text workflow lines (PARSEFRAGPIPEWORKFLOW msbooster.config)
 * @return numThreads as an Integer
 */
def msboosterThreads(params_text) {
    def lines = params_text.toString().readLines().findAll { line -> line ==~ /numThreads = \d+/ }
    if (lines.size() != 1) {
        error("msboosterThreads: the MSBooster params need exactly one numThreads line (workflow.threads or --fragpipe_threads), got ${lines}")
    }
    return lines[0].substring('numThreads = '.length()) as Integer
}

/**
 * The PIN paths MSBooster is given, relative to FragPipe's work directory: <group>/<name>, or the bare name
 * for the blank group (whose files sit in the work directory, CmdMSBooster.java:203-219).
 *
 * @param pin_names staged PIN file names
 * @param pin_dirs  FragPipe group of each PIN, in the same order ('' for the blank group)
 * @return PIN paths in the same order
 */
def msboosterPinPaths(pin_names, pin_dirs) {
    if (pin_names.size() != pin_dirs.size()) {
        error("msboosterPinPaths: ${pin_names.size()} PINs but ${pin_dirs.size()} groups")
    }
    if (pin_dirs.any { dir -> dir == null }) {
        error("msboosterPinPaths: a PIN has no FragPipe group (meta.group)")
    }
    return [pin_names, pin_dirs].transpose().collect { name, dir -> dir ? "${dir}/${name}".toString() : name.toString() }
}

/**
 * MSBooster's output directory, relative to the task directory: MSBooster under the deepest directory common to
 * every PIN's parent (one PIN: its parent; Constants.updateOutputDirectory, MyFileUtils.findDeepestCommonDirectory).
 * It holds spectraRT.tsv, spectraRT_full.tsv, spectraRT.predicted.bin and MSBooster_plots/.
 *
 * @param pin_paths msboosterPinPaths() of the run
 * @return 'MSBooster' or '<dir>/MSBooster'
 */
def msboosterOutputDir(pin_paths) {
    def parents = pin_paths.collect { path -> path.contains('/') ? path.substring(0, path.lastIndexOf('/')).tokenize('/') : [] }
    def common = parents.inject { a, b -> [a, b].transpose().takeWhile { x, y -> x == y }.collect { x, _y -> x } }
    return (common + ['MSBooster']).join('/')
}

/**
 * Shell lines that replace each staged PIN link by a copy and move grouped PINs to their group path; the
 * blank group's PINs stay put. MSBooster writes <pin>_edited.pin next to the PIN's canonical path
 * (PercolatorFormatter.java:190, 292-300), so the PIN must be a file in the task dir, not a link.
 *
 * @param pin_paths msboosterPinPaths() of the run
 * @return shell lines
 */
def msboosterPinStaging(pin_paths) {
    def copies = linkCopies(pin_paths.collect { path -> path.substring(path.lastIndexOf('/') + 1) })
    return copies + '\n' + pin_paths.findAll { path -> path.contains('/') }.collect { path ->
        "mkdir -p ${shellArgs([path.substring(0, path.lastIndexOf('/'))])} && mv ${shellArgs([path.substring(path.lastIndexOf('/') + 1)])} ${shellArgs([path])}"
    }.join('\n')
}

/**
 * msbooster_params.txt as CmdMSBooster.configure writes it (CmdMSBooster.java:160-212): the workflow
 * lines of parse_fragpipe_workflow.py without useIM, findBestImModel and numThreads, then the first two ANDed with
 * "the run has timsTOF data" (CmdMSBooster.java:151-157, 170, 178), then the runtime keys numThreads (the
 * parser's, which must equal threads; MSBooster's output depends on it),
 * DiaNN (FragPipe's bundled DIA-NN, Diann.fallbackDiannPath), fragger, unimodObo, mzmlDirectory and
 * pinPepXMLDirectory. The keys are written in CmdMSBooster's order. Paths in the task directory use {TASK_DIR}.
 *
 * @param params_text    workflow lines (PARSEFRAGPIPEWORKFLOW msbooster.config)
 * @param tims_tof       true when an LC-MS file of the run is timsTOF data
 * @param threads        FragPipe's thread count (workflow.threads / --threads), the parser's numThreads
 * @param tools          fragpipeTools() map
 * @param fragger_params staged MSFragger params file name MSBooster reads
 * @param pins           PIN paths relative to the task directory, in FragPipe's order (msboosterPinPaths)
 * @param library        the staged msbooster.spectral-library-path file (spectral_library/<name>), [] without one
 * @return msbooster_params.txt text
 */
def msboosterParamsText(params_text, tims_tof, threads, Map tools, fragger_params, pins, library) {
    def library_files = [library].flatten().findAll { f -> f }
    def has_token = params_text.toString().contains('{SPECTRAL_LIBRARY}')
    // CmdMSBooster.java:180-184: the three prediction files name the library only when libraryPath is not empty
    if (has_token != !library_files.isEmpty()) {
        error("msboosterParamsText: the workflow params ${has_token ? 'name' : 'do not name'} a spectral library but ${library_files.size()} was staged")
    }
    def lines = params_text.toString().replace('{SPECTRAL_LIBRARY}', library_files ? "{TASK_DIR}/${library_files[0]}".toString() : '').readLines()
    // numThreads is the parser's (workflow.threads / --threads): MSBooster's RT loess output depends on it,
    // and the task is sized from the same count (conf/modules.config)
    def thread_lines = lines.findAll { line -> line.startsWith('numThreads = ') }
    if (thread_lines != ["numThreads = ${threads}".toString()]) {
        error("msboosterParamsText: the MSBooster workflow params must hold numThreads = ${threads}, got ${thread_lines}")
    }
    def im_keys = ['useIM', 'findBestImModel', 'numThreads']
    def kept = lines.findAll { line -> !im_keys.any { k -> line.startsWith("${k} = ") } }
    def im_lines = im_keys.take(2).collect { k ->
        def line = lines.find { l -> l.startsWith("${k} = ") }
        if (line == null) {
            error("msboosterParamsText: ${k} missing from the MSBooster workflow params")
        }
        "${k} = ${tims_tof && line.substring(k.length() + 3) == 'true' ? 'true' : 'false'}"
    }
    def pin_names = (pins instanceof List ? pins : [pins]).collect { p -> p.toString() }
    // CmdMSBooster.java:200-212: an empty HashSet filled by add() in LC-MS order. new HashSet(Collection) presizes
    // the table (32 buckets for 12 PINs on JDK 17, 16 in FragPipe's set), which changes the iteration order.
    def pin_set = new HashSet<String>()
    pin_names.each { name -> pin_set.add(name) }
    def runtime = [
        "numThreads = ${threads}",
        "DiaNN = ${tools.diann}",
        "fragger = {TASK_DIR}/${fragger_params}",
        "unimodObo = ${tools.unimod_obo}",
        // CmdMSBooster.java:196-199: each LC-MS directory is followed by a space
        'mzmlDirectory = {TASK_DIR} ',
        // CmdMSBooster.java:203-219: a HashSet of the PIN paths, written in its iteration order
        "pinPepXMLDirectory = ${pin_set.join(' ')}",
    ]
    // FragPipe's write order (CmdMSBooster.java:160-212); spectraPredFile/RTPredFile/IMPredFile only with a library
    def key_order = ['useDetect', 'numThreads', 'DiaNN', 'renamePin', 'useRT', 'useSpectra', 'useIM', 'fragger', 'deletePreds',
        'rtModel', 'spectraModel', 'imModel', 'findBestRtModel', 'findBestSpectraModel', 'findBestImModel', 'KoinaURL',
        'spectraPredFile', 'RTPredFile', 'IMPredFile', 'rtSearchModelsString', 'ms2SearchModelsString', 'imSearchModelsString',
        'FragmentationType', 'unimodObo', 'mzmlDirectory', 'pinPepXMLDirectory']
    def by_key = (kept + im_lines + runtime).collectEntries { line -> [(line.substring(0, line.indexOf(' = '))): line] }
    def unknown = by_key.keySet() - key_order
    if (unknown) {
        error("msboosterParamsText: keys without a position in CmdMSBooster's write order: ${unknown}")
    }
    return key_order.findAll { k -> by_key.containsKey(k) }.collect { k -> "${by_key[k]}\n" }.join('')
}

/**
 * The files and options of an IonQuant run as CmdIonquant.configure writes them (CmdIonquant.java:178-337),
 * with the task directory ({TASK_DIR}) as FragPipe's work dir.
 *
 * - groups: FragPipe's blank group "" first (its files sit in the work-dir root, LcmsFileGroup.outputDir),
 *   then the named groups sorted by name like FragPipe's TreeMap
 * - filelist: filelist_ionquant.txt, "flag\tvalue", one --psm per group relative to the work dir, then
 *   --specdir with the absolute LC-MS directory (CmdIonquant.java:310-337)
 * - annotations: for the isobaric pass (--perform-isoquant 1), one "<channel> <sample_name>" file per group
 *   from the multi-plex TMT annotation, <group>/<group>_annotation.txt (_annotation.txt for the blank
 *   group; TmtiPanel.java:712-713), passed as --annotation <psm>=<annotation> (CmdIonquant.java:204-219)
 * - args: the pass's options with the --annotation pairs before --minexps, then --multidir . when a
 *   group is named (CmdIonquant.java:221-222, 324-327)
 *
 * @param group_names      staged group directory names
 * @param root_group_dir   the staged directory of the blank group ('' when there is none)
 * @param config_cli       IonQuant options of the pass from the workflow parser, in CmdIonquant order
 * @param annotation_lines lines of the multi-plex TMT annotation (plex, channel, sample, sample_name, ...), [] if none
 * @param spectra_dir      directory of the staged LC-MS files, under {TASK_DIR}
 * @return Map with named_groups, filelist (text), annotations ([path, text] pairs) and args (List)
 */
def ionquantRun(group_names, root_group_dir, config_cli, annotation_lines, spectra_dir) {
    def root = (root_group_dir ?: '').toString()
    def named = group_names.collect { n -> n.toString() }.findAll { n -> n != root }.sort()
    if (!root && !named) {
        error("IonQuant: no psm.tsv group directory was staged")
    }
    def options = config_cli.toString().tokenize()
    def option_value = { flag -> def i = options.indexOf(flag); i >= 0 && i + 1 < options.size() ? options[i + 1] : null }
    // CmdIonquant.java:288-294: MBR needs named experiments
    if (!named && option_value.call('--mbr') == '1') {
        error("IonQuant with MBR requires designating LCMS runs to experiments. Please assign the experiment groups in the 'Workflow' tab, and make sure that the 'Spec Lib' tab is disabled.")
    }
    def plexes = (root ? [''] : []) + named
    def psm = { plex -> plex ? "${plex}/psm.tsv".toString() : 'psm.tsv' }
    def filelist = (['flag\tvalue'] + plexes.collect { plex -> "--psm\t${psm.call(plex)}" } + ["--specdir\t${spectra_dir}"])
        .collect { line -> "${line}\n" }.join('')
    def annotations = []
    def annotation_args = []
    if (option_value.call('--perform-isoquant') == '1') {
        def rows = (annotation_lines ?: []).drop(1).collect { line -> line.toString().split('\t', -1) as List }
        if (!rows) {
            error("--perform-isoquant 1 needs the TMT annotation (TMT-Integrator annotations, FragpipeRun.java:2005-2031)")
        }
        plexes.each { plex ->
            def lines = rows.findAll { cols -> cols[0] == plex }.collect { cols -> "${cols[1]} ${cols[3]}\n" }
            if (!lines) {
                error("TMT annotation has no rows for group '${plex}'")
            }
            def path = plex ? "${plex}/${plex}_annotation.txt".toString() : '_annotation.txt'
            if ((psm.call(plex) + path).contains('=')) {
                error("There are '=' sign in ${psm.call(plex)} or ${path}. Please rename the file or directory.")
            }
            annotations << [path, lines.join('')]
            annotation_args += ['--annotation', "{TASK_DIR}/${psm.call(plex)}={TASK_DIR}/${path}".toString()]
        }
    }
    def minexps = options.indexOf('--minexps')
    if (minexps < 0) {
        error("IonQuant options lack --minexps (CmdIonquant.java:216-219 always sets it): ${config_cli}")
    }
    def args = options.take(minexps) + annotation_args + options.drop(minexps) + (named ? ['--multidir', '.'] : [])
    return [named_groups: named, filelist: filelist, annotations: annotations, args: args]
}

/**
 * Shell lines that lay out an IonQuant task directory as FragPipe's work dir: the blank group's files in
 * the root (LcmsFileGroup.outputDir), named groups as real directories (IonQuant updates them in place).
 *
 * @param root_group_dir staged blank-group directory ('' when there is none)
 * @param named_groups   staged named group directories
 * @return shell lines
 */
def ionquantStaging(root_group_dir, named_groups) {
    def lines = root_group_dir ? ["cp -rL '${root_group_dir}/.' ."] : []
    lines += named_groups.collect { g -> "if [ -L '${g}' ]; then cp -rL '${g}' '${g}.copy' && rm '${g}' && mv '${g}.copy' '${g}'; fi" }
    return lines.collect { line -> "    ${line}" }.join('\n')
}

/**
 * Shell lines that write the files of an IonQuant run (ionquantRun): filelist_ionquant.txt, the per-group
 * TMT annotations, modmasses_ionquant.txt (one mass per line, written even when empty; CmdIonquant.java:340-350)
 * and experiment_annotation.tsv (CmdIonquant.java:366-377).
 *
 * @param run                   ionquantRun() result
 * @param modmasses             mod mass lines, or null for no --modlist
 * @param experiment_annotation experiment_annotation.tsv lines ([] = none)
 * @return shell lines
 */
def ionquantFiles(Map run, modmasses, experiment_annotation) {
    def lines = [writeTextFile('filelist_ionquant.txt', run.filelist)]
    lines += run.annotations.collect { a -> writeTextFile(a[0], a[1]) }
    if (modmasses != null) {
        lines << writeTextFile('modmasses_ionquant.txt', modmasses.collect { m -> "${m}\n" }.join(''))
    }
    if (experiment_annotation) {
        lines << writeTextFile('experiment_annotation.tsv', experiment_annotation.collect { l -> "${l}\n" }.join(''))
    }
    return lines.join('\n')
}

/**
 * Shell lines that copy an IonQuant work dir into <prefix>/ in FragPipe's layout: the root's tables and
 * the files IonQuant read (not experiment_annotation.tsv, the module's own output), the named groups; and
 * the blank group's updated files (plus ion.tsv) into ionquant_root/<root_group_dir>/ for a later pass.
 *
 * @param prefix         task output directory
 * @param root_group_dir staged blank-group directory ('' when there is none)
 * @param named_groups   named group directories
 * @return shell lines
 */
def ionquantCollect(prefix, root_group_dir, named_groups) {
    def lines = [
        'shopt -s nullglob',
        "for f in *.tsv *.csv filelist_ionquant.txt modmasses_ionquant.txt *_annotation.txt; do [ \"\$f\" = experiment_annotation.tsv ] || [ ! -e \"\$f\" ] || cp -L \"\$f\" '${prefix}/'; done",
        'shopt -u nullglob',
    ]
    if (named_groups) {
        lines << "cp -rL ${named_groups.collect { g -> "'${g}'" }.join(' ')} '${prefix}/'"
    }
    if (root_group_dir) {
        lines << "mkdir -p 'ionquant_root/${root_group_dir}'"
        lines << "for f in '${root_group_dir}'/*; do cp -L \"\$(basename \"\$f\")\" 'ionquant_root/${root_group_dir}/'; done"
        lines << "[ ! -f ion.tsv ] || cp -L ion.tsv 'ionquant_root/${root_group_dir}/'"
    }
    return lines.collect { line -> "    ${line}" }.join('\n')
}

/**
 * FP-Meta DbOptimizer argv after the jar (CmdMetaproteomics.java:70-111): FragPipe passes its work dir as
 * both --projectDir and --outdir, cmd-line-opts as one last element when not blank.
 *
 * @param values      parse_fragpipe_workflow.py metaproteomics values (qvalue, deltaHyperscore, ...)
 * @param project_dir FragPipe work dir for FP-Meta
 * @param fasta       database path
 * @param decoy_tag   decoy prefix (FragpipeRun.java:1442)
 * @param names_dmp   NCBI names.dmp path
 * @param nodes_dmp   NCBI nodes.dmp path
 * @return argv words
 */
def metaproteomicsArgv(Map values, project_dir, fasta, decoy_tag, names_dmp, nodes_dmp) {
    def keys = ['qvalue', 'deltaHyperscore', 'minPeptCntPerProt', 'minUniqPeptCntPerProt', 'minUniqPeptCnt', 'hostName', 'iterations', 'cmdLineOpts']
    def missing = keys.findAll { k -> values[k] == null }
    if (missing) {
        error("metaproteomicsArgv: parse_fragpipe_workflow.py metaproteomics values lack ${missing}")
    }
    def argv = ['DbOptimizer', '--projectDir', project_dir, '--fastaFile', fasta, '--outdir', project_dir, '--decoyTag', decoy_tag]
    argv += ['qvalue', 'deltaHyperscore', 'minPeptCntPerProt', 'minUniqPeptCntPerProt', 'minUniqPeptCnt'].collectMany { k -> ["--${k}", values[k]] }
    argv += ['--taxonNameFile', names_dmp, '--taxonNodeFile', nodes_dmp, '--hostName', values.hostName, '--iterations', values.iterations]
    if (values.cmdLineOpts) {
        argv << values.cmdLineOpts
    }
    return argv.collect { w -> w.toString() }
}

/**
 * Shell lines that copy the tables FPOP coADAPTr edits into <outdir> in FragPipe's layout (symlinks
 * dereferenced, so the in-place edits of FPOP_PSMwriter / FPOP_DiannReport_Editor leave the upstream
 * files untouched; CmdFPOPcoadaptr.java:29-58, 429-433).
 *
 * Layouts: 'tree' copies each directory whole (IonQuant results: every psm.tsv inside); 'groups' copies
 * <name>/psm.tsv (Philosopher per-group results); 'diann' copies <name>/report.tsv (DIA-NN dia-quant-output).
 *
 * @param layout tree, groups or diann
 * @param staged staged input paths, e.g. input/1/<name>
 * @param outdir FragPipe work dir for the editor
 * @return shell lines
 */
def fpopCoadaptrStaging(layout, staged, outdir) {
    def tables = [tree: null, groups: 'psm.tsv', diann: 'report.tsv']
    if (!tables.containsKey(layout)) {
        error("FPOP coADAPTr: unknown layout '${layout}' (tree, groups, diann)")
    }
    def paths = (staged instanceof List ? staged : [staged]).collect { p -> p.toString() }
    if (!paths) {
        error("FPOP coADAPTr: nothing staged")
    }
    def names = paths.collect { p -> p.tokenize('/').last() }
    def dups = names.findAll { n -> names.count(n) > 1 }.unique()
    if (dups) {
        error("FPOP coADAPTr: two inputs are named ${dups}")
    }
    def table = tables[layout]
    def lines = ["mkdir -p '${outdir}'"]
    [paths, names].transpose().each { p, n ->
        if (table == null) {
            lines << "cp -rL '${p}' '${outdir}/${n}'"
        } else {
            lines << "mkdir -p '${outdir}/${n}'"
            lines << "cp -L '${p}/${table}' '${outdir}/${n}/${table}'"
        }
    }
    return lines.join('\n')
}

/**
 * The psm.tsv files FragAnnotator's export rewrites, found from fragpipe-files.fp-manifest as
 * ResultProcessor.processManifestFile/goThroughFolder does: only the work-dir root when any experiment
 * is blank or SpecLibGen runs, else <experiment>/ or <experiment>_<bioreplicate>/ of every line.
 * The EXPORTMATCHEDFRAGMENTS stub edits these files in place of the export.
 *
 * @param manifest_text fragpipe-files.fp-manifest text
 * @param speclibgen    true when SpecLibGen runs (speclibgen.run-speclibgen)
 * @return paths under fragpipe/, sorted
 */
def matchedFragmentsExportPsms(manifest_text, speclibgen) {
    def rows = manifest_text.toString().readLines().findAll { line -> line }.collect { line -> line.split('\t', -1) as List }
    if (speclibgen || rows.any { cols -> cols[1] == '' }) {
        return ['fragpipe/psm.tsv']
    }
    return rows.collect { cols -> cols[2] == '' ? cols[1] : "${cols[1]}_${cols[2]}" }.unique().sort().collect { d -> "fragpipe/${d}/psm.tsv".toString() }
}

/**
 * Shell lines that build the FragPipe work dir FragmentsExportMainClass reads (FragAnnotator-1.0.0.jar,
 * ResultProcessor + ExportFragments) in ./fragpipe/: the manifest, fragpipe.workflow with
 * speclibgen.run-speclibgen (the only workflow key the export reads, ResultProcessor.processWorkflow),
 * the named groups as directories and the blank group's tables in the root (inner_defined_empty_exp),
 * and umpire-se.params when DIA-Umpire runs. Group directories are copied: the export rewrites psm.tsv.
 *
 * @param manifest       staged manifest name
 * @param group_names    staged group directory names (under groups/)
 * @param root_group_dir the blank group's stand-in name ('' when there is none)
 * @param diaumpire      true when DIA-Umpire runs
 * @param speclibgen     true when SpecLibGen runs
 * @return shell lines
 */
def matchedFragmentsStaging(manifest, group_names, root_group_dir, diaumpire, speclibgen) {
    def lines = ['mkdir -p fragpipe', "cp '${manifest}' fragpipe/fragpipe-files.fp-manifest",
        writeTextFile('fragpipe/fragpipe.workflow', "speclibgen.run-speclibgen=${speclibgen ? 'true' : 'false'}\n")]
    group_names.each { n ->
        lines << (n == root_group_dir ? "cp -rL 'groups/${n}/.' fragpipe/" : "cp -rL 'groups/${n}' 'fragpipe/${n}'")
    }
    if (diaumpire) {
        lines << 'touch fragpipe/umpire-se.params'
    }
    return lines.join('\n')
}

/**
 * Shell lines that copy every group's tables from ./fragpipe/ to results/<group dir>/ (the blank group
 * from its staged directory plus the root psm.tsv). A group the export did not find keeps its psm.tsv.
 *
 * @param group_names    staged group directory names
 * @param root_group_dir the blank group's stand-in name ('' when there is none)
 * @return shell lines
 */
def matchedFragmentsCollect(group_names, root_group_dir) {
    def lines = ['mkdir -p results']
    group_names.each { n ->
        if (n == root_group_dir) {
            lines << "cp -rL 'groups/${n}' 'results/${n}'"
            lines << "cp fragpipe/psm.tsv 'results/${n}/psm.tsv'"
        } else {
            lines << "cp -rL 'fragpipe/${n}' 'results/${n}'"
        }
    }
    return lines.join('\n')
}

/**
 * DIA-NN inputs largest first, as CmdDiann orders filelist_diann.txt (CmdDiann.java:268-294, 389-397): a .d
 * directory counts every regular file under it. Equal sizes keep the given (LC-MS table) order; FragPipe keeps
 * them in one TreeMap bucket in Set<Path> iteration order, which is not reproducible. Reads file sizes.
 *
 * @param files DIA-NN input files and .d directories
 * @return the files, largest first
 */
def diannSizeOrder(files) {
    def size = { java.nio.file.Path p ->
        if (!java.nio.file.Files.isDirectory(p)) {
            return java.nio.file.Files.size(p)
        }
        def total = 0L
        java.nio.file.Files.walk(p, java.nio.file.FileVisitOption.FOLLOW_LINKS).withCloseable { stream ->
            total = stream.filter { f -> java.nio.file.Files.isRegularFile(f) }.mapToLong { f -> java.nio.file.Files.size(f) }.sum()
        }
        return total
    }
    def list = files instanceof List ? files : [files]
    return list.withIndex().collect { f, i -> [f, size.call(f), i] }
        .sort { a, b -> b[1] <=> a[1] ?: a[2] <=> b[2] }
        .collect { item -> item[0] }
}

/**
 * DIA-NN argv of CmdDiann.java:313-409 (DIA-NN 1.8.2 Beta 8 branches), cwd = the group work dir, split around
 * the place of task.ext.args: head (binary to protein-inference and args_post), tail (--cfg; its "--" suffix
 * is FragPipe's on Linux).
 *
 * @param cfg      parsed diann entry (args, args_post, cmd_opts, library, run_plex, redo_protein_inference)
 * @param diann    DIA-NN binary
 * @param threads  thread count (the task's CPUs)
 * @param library  staged external library path (cfg.library set), e.g. library_input/x.tsv
 * @param fasta    staged database name
 * @return Map with head and tail argv words
 */
def diannArgv(Map cfg, diann, threads, library, fasta) {
    // --lib (CmdDiann.java:315-326): plexDIA library, external library, else library.tsv
    def lib = cfg.run_plex ? 'library_2.tsv' : (cfg.library ? "{TASK_DIR}/${library}" : 'library.tsv')
    // Protein inference for DIA-NN < 2.0 (CmdDiann.java:345-352)
    def prot_inf = cfg.redo_protein_inference ? ['--relaxed-prot-inf', '--fasta', "{TASK_DIR}/${fasta}"] : ['--no-prot-inf']
    def head = [diann, '--lib', lib, '--threads', threads, '--verbose', '1', '--out', 'dia-quant-output/report.tsv'] +
        (cfg.args ?: '').toString().tokenize(' ') + prot_inf + (cfg.args_post ?: '').toString().tokenize(' ') +
        // diann.cmd-opts is one argv element (cmd.add(additionalCmdOpts), CmdDiann.java:386-388)
        (cfg.cmd_opts ? [cfg.cmd_opts] : [])
    return [head: head.collect { w -> w.toString() }, tail: ['--cfg', '{TASK_DIR}/filelist_diann.txt--']]
}

/**
 * filelist_diann.txt: "--f <path>" per DIA-NN input in the given order (CmdDiann.java:389-397).
 *
 * @param spectra staged DIA-NN inputs (<name> in the task dir), largest first (diannSizeOrder)
 * @return file text
 */
def diannFilelistText(spectra) {
    def list = spectra instanceof List ? spectra : [spectra]
    if (!list) {
        error("DIA-NN: no input files")
    }
    return list.collect { f -> "--f {DIANN_INPUT_DIR}/${f.toString().tokenize('/').last()}\n" }.join('')
}

/**
 * Shell lines that set DIANN_INPUT_DIR, the absolute directory DIA-NN reads its staged inputs from.
 *
 * DIA-NN 1.8.2 Beta 8 keeps a run's path in a 256-byte buffer: a path of 256 characters or more aborts it
 * with "buffer overflow detected" (measured with FragPipe 24.0's binary: 255 runs, 256 aborts; long --lib
 * and --cfg paths are fine). FragPipe passes absolute paths (CmdDiann.java:389-397), so the task directory
 * is used as it is when every input path fits, which keeps filelist_diann.txt and the reports as
 * FragPipe writes them. For a longer work dir (S3/Fusion, deep scratch) the inputs are read through a
 * symlink in a directory of the task's $TMPDIR, else /tmp (mktemp -d, so tasks sharing a host never collide), removed
 * when the task script exits. DIA-NN and FragReporter then record that directory as the runs' location
 * (report File.Name, matrix and abundance run columns, log), in place of the task dir: the same
 * results with a different LC-MS directory prefix, as with any other work dir.
 *
 * @param spectra staged DIA-NN inputs (<name> in the task dir)
 * @return shell lines
 */
def diannInputDirLines(spectra) {
    def list = spectra instanceof List ? spectra : [spectra]
    def longest = list.collect { f -> f.toString().tokenize('/').last().length() }.max()
    def max_path = 255
    return [
        'DIANN_INPUT_DIR="$PWD"',
        "if [ \$(( \${#DIANN_INPUT_DIR} + 1 + ${longest} )) -gt ${max_path} ]; then",
        '    DIANN_LINK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/diann.XXXXXX")"',
        '    trap \'rm -rf "$DIANN_LINK_DIR"\' EXIT',
        '    DIANN_INPUT_DIR="$DIANN_LINK_DIR/s"',
        '    ln -s "$PWD" "$DIANN_INPUT_DIR"',
        'fi',
    ].join('\n')
}

/**
 * The FragPipe headless work layout: LC-MS file moves, manifest.fp-manifest and TMT annotations.
 *
 * manifest.fp-manifest has one line per LC-MS file (path, experiment, bioreplicate, data type), the four
 * columns FragPipe's manifest loader reads (TabWorkflow.java:1255-1282). With a TMT annotation FragPipe
 * looks for exactly one *annotation.txt next to the LC-MS files of a group (TmtiPanel.java:868-880), so each
 * named group's files move into raw_files/<group>/ with <group>_annotation.txt ("<channel> <sample>" lines);
 * the blank group "" stays in raw_files/ with _annotation.txt, FragPipe's default name (TmtiPanel.java:712-720).
 *
 * @param rows            [file name, experiment, bioreplicate, data type, group] per LC-MS file
 * @param annotation_rows [group, channel, sample name] per TMT channel ([] for label-free)
 * @return Map with moves ([from, to] pairs), manifest (text, paths under {TASK_DIR}) and annotations ([path, text] pairs)
 */
def headlessLayout(rows, annotation_rows) {
    def is_tmt = !annotation_rows.isEmpty()
    def moves = []
    def manifest = rows.collect { row ->
        def name = row[0]
        def group = row[4]
        def path = "raw_files/${name}"
        if (is_tmt) {
            if (group == null) {
                error("FragPipe headless: no FragPipe group for ${name}")
            }
            if (group) {
                moves << [path.toString(), "raw_files/${group}/${name}".toString()]
                path = "raw_files/${group}/${name}"
            }
        }
        "{TASK_DIR}/${path}\t${row[1] ?: ''}\t${row[2] == null ? '' : row[2]}\t${row[3]}\n"
    }.join('')
    def annotations = annotation_rows.collect { a -> a[0].toString() }.unique().sort().collect { group ->
        def dir = group ? "raw_files/${group}" : 'raw_files'
        def text = annotation_rows.findAll { a -> a[0].toString() == group }.collect { a -> "${a[1]} ${a[2]}\n" }.join('')
        ["${dir}/${group}_annotation.txt".toString(), text]
    }
    return [moves: moves, manifest: manifest, annotations: annotations]
}

/**
 * Stops a FragPipe headless run whose manifest rows do not describe exactly the staged LC-MS files
 * (nf-core/msproteomics a3f8742): rows are matched to the files by exact name, never by substring, and a
 * staged file without a row, a row naming no staged file, a file listed twice, no rows at all or a TMT
 * group without annotation channels are errors, not silently dropped files.
 *
 * @param raw_files       staged LC-MS files (one path or a list)
 * @param rows            [file name, experiment, bioreplicate, data type, group] per LC-MS file
 * @param annotation_rows [group, channel, sample name] per TMT channel ([] for label-free)
 * @return true when the rows match the files
 */
def checkHeadlessRows(raw_files, rows, annotation_rows) {
    def staged = (raw_files instanceof List ? raw_files : [raw_files]).collect { f -> f.toString().tokenize('/').last() }
    if (rows.isEmpty()) {
        error("FragPipe headless: no manifest rows for the ${staged.size()} staged LC-MS file(s)")
    }
    def names = rows.collect { row -> row[0].toString() }
    def duplicates = names.findAll { n -> names.count(n) > 1 }.unique()
    if (duplicates) {
        error("FragPipe headless: LC-MS file(s) listed more than once in the manifest rows: ${duplicates.join(', ')}")
    }
    def without_row = staged - names
    if (without_row) {
        error("FragPipe headless: staged LC-MS file(s) without a manifest row: ${without_row.join(', ')}")
    }
    def unstaged = names - staged
    if (unstaged) {
        error("FragPipe headless: manifest row(s) naming no staged LC-MS file: ${unstaged.join(', ')}")
    }
    if (!annotation_rows.isEmpty()) {
        def annotated = annotation_rows.collect { a -> a[0].toString() }.unique()
        def missing = rows.findAll { row -> row[4] != null && !(row[4].toString() in annotated) }.collect { row -> "${row[0]} (${row[4]})" }
        if (missing) {
            error("FragPipe headless: TMT file(s) whose group has no annotation channels: ${missing.join(', ')}")
        }
    }
    return true
}

/**
 * Shell lines that lay out the FragPipe headless task directory (headlessLayout): staged .raw links
 * replaced by copies (rawLinkCopies), TMT group moves, manifest.fp-manifest and the annotations.
 *
 * @param raw_files staged LC-MS files (raw_files/<name>)
 * @param layout    headlessLayout() result
 * @return shell lines
 */
def headlessStagingLines(raw_files, Map layout) {
    def lines = [rawLinkCopies(raw_files)] + layout.moves.collect { m -> "mkdir -p '${m[1].substring(0, m[1].lastIndexOf('/'))}' && mv '${m[0]}' '${m[1]}'" }
    lines << writeTextFile('manifest.fp-manifest', layout.manifest)
    lines += layout.annotations.collect { a -> writeTextFile(a[0], a[1]) }
    return lines.findAll { l -> l }.join('\n')
}

/**
 * The .workflow FragPipe headless runs (run.workflow): the user's workflow with database.db-path set to the
 * staged database (replaced, or appended when absent), which FragPipe has no flag for. The tools come
 * from FragPipe's own flags and lookups, as in a user's headless run: --config-tools-folder (MSFragger,
 * IonQuant, diaTracer), FragPipeMain.getPhilosopherBin and the built-in DIA-NN (Diann.java:57); no
 * component reads a tool path key from the workflow, and extra keys would change the HashMap FragPipe
 * loads the workflow from, which decides both-true radio twins (parse_fragpipe_workflow.py _resolve_radio_twins).
 *
 * @param workflow_text the user's .workflow text
 * @param database      staged database name
 * @return run.workflow text
 */
def headlessWorkflowText(workflow_text, database) {
    def db_line = "database.db-path={TASK_DIR}/${database}".toString()
    def lines = workflow_text.toString().readLines()
    lines = lines.any { l -> l.startsWith('database.db-path=') } ? lines.collect { l -> l.startsWith('database.db-path=') ? db_line : l } : lines + [db_line]
    return lines.collect { l -> "${l}\n" }.join('')
}

/**
 * Java heap of a FragPipe JVM tool in GB: ext.java_xmx, else the task memory. FragPipe has no fixed
 * default (its RAM setting, else OsUtils.getDefaultXmx(), FragpipeRun.java:1308), so a task without
 * either stops.
 *
 * @param task Task (ext.java_xmx, memory)
 * @return heap size in GB
 */
def fragpipeJavaMemGb(task) {
    if (task.ext.java_xmx) {
        return task.ext.java_xmx as int
    }
    if (!task.memory) {
        error("${task.process}: neither ext.java_xmx nor task.memory is set; the Java heap (-Xmx) is the task memory, as FragPipe's RAM setting (FragpipeRun.java:1308)")
    }
    return task.memory.toGiga()
}

/**
 * TMT label-check thresholds of tmt_qc.py: ext.pass_threshold (0.95), ext.warn_threshold (0.85) and
 * ext.prob_threshold (0.95).
 *
 * @param task Task (ext.pass_threshold, ext.warn_threshold, ext.prob_threshold)
 * @return Map with pass, warn and prob
 */
def tmtLabelcheckThresholds(task) {
    return [
        pass: task.ext.pass_threshold ?: 0.95,
        warn: task.ext.warn_threshold ?: 0.85,
        prob: task.ext.prob_threshold ?: 0.95
    ]
}

/**
 * Whether Philosopher Report is skipped after Filter (ext.skip_report, FreeQuant / LabelQuant update the
 * workspace before Report: FragpipeRun.java:2463-2471).
 *
 * @param task Task (ext.skip_report)
 * @return true when Report is skipped
 */
def philosopherSkipReport(task) {
    return task.ext.skip_report ?: false
}

/**
 * Names of the TMT-Integrator work dir (ext.work_dir_name, default the prefix), report dir
 * (ext.report_dir_name, default tmt-report) and config file (ext.config_file_name, default
 * tmt-integrator-conf.yml). FragPipe's second run for FPOP TMT writes tmt-report-unmod with
 * tmt-integrator-conf-unmod.yml (CmdTmtIntegrator.java:62, :149-152).
 *
 * @param task Task (ext.work_dir_name, ext.report_dir_name, ext.config_file_name)
 * @param prefix default work dir name
 * @return Map with work_dir, report_dir and config_name
 */
def tmtIntegratorNames(task, prefix) {
    return [
        work_dir   : task.ext.work_dir_name ?: prefix,
        report_dir : task.ext.report_dir_name ?: 'tmt-report',
        config_name: task.ext.config_file_name ?: 'tmt-integrator-conf.yml'
    ]
}

/**
 * Name of a Crystal-C params file: <ccParamsFilePrefix>-<index>-<pepXML name>.params. The index is counted
 * over the whole run by fragpipe_search (meta.crystalc_index, CmdCrystalc.java:176-185).
 *
 * @param meta Sample meta map (crystalc_index)
 * @param pepxml_file pepXML file the run reads
 * @return params file name
 */
def crystalcParamsName(meta, pepxml_file) {
    if (meta.crystalc_index == null) {
        error("CRYSTALC needs meta.crystalc_index (the run-wide index fragpipe_search assigns, CmdCrystalc.java:176-185)")
    }
    return "crystalc-${meta.crystalc_index}-${pepxml_file.name}.params".toString()
}

/**
 * MSFragger --partial argument of one chunk of a split-database search (meta.chunk_id).
 *
 * @param meta Sample meta map (chunk_id)
 * @return ['--partial', chunk_id], or [] without a chunk
 */
def msfraggerPartialArgs(meta) {
    return meta.chunk_id != null ? ['--partial', "${meta.chunk_id}".toString()] : []
}

/**
 * Name Philosopher LabelQuant requires for the TMT annotation of a group (CmdLabelquant.java:86-91).
 *
 * @param meta Group meta map (group)
 * @return <group>_annotation.txt
 */
def labelquantAnnotationName(meta) {
    return "${meta.group}_annotation.txt".toString()
}

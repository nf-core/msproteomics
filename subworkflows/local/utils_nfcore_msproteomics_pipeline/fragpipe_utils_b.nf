/*
 * Command builders for the FragPipe validation, inference, quantification and export modules:
 * each returns the exact argv (and derived names) FragPipe 24.0 runs, so the modules only run them.
 * Pure functions, no I/O; unit tests in tests/fragpipe_utils_b.nf.test.
 */

include { expandJarTokens; fragpipeTools } from './fragpipe_utils_a'

/**
 * Single-quote one token for bash; a fragpipeJarTokens token (a jar the task finds) becomes its
 * double-quoted shell variable inside the word.
 *
 * @param token Value to quote
 * @return The token as one bash word
 */
def shellQuote(token) {
    return expandJarTokens("'" + token.toString().replace("'", "'\\''") + "'", "'\"", "\"'")
}

/**
 * Percolator and PercolatorOutputToPepXML argv exactly as CmdPercolator builds them.
 *
 * Basename: the pepXML stem without _rank<N> (CmdPercolator.java:133-137); with Crystal-C the
 * PIN and TSV names drop the trailing _c (CmdPercolator.java:145-150). Percolator: cmd-opts,
 * --num-threads, results flags (--only-psms picks PSM or peptide flags), decoy pattern, then
 * --weights <stripped>_percolator.weights with Write sub mzML (CmdPercolator.java:180-183, writeSubMzml =
 * tabRun.isWriteSubMzml(), FragpipeRun.java:1642), --init-weights <weights> --static --override for a second-pass
 * _sub.mzML input (CmdPercolator.java:185-201; {TASK_DIR}/weights/<name>, see shellArgs), then
 * the edited PIN when MSBooster ran, else the MSFragger PIN (CmdPercolator.java:152-213). The
 * binary is tools/percolator_3_7_1/linux/percolator (CmdPercolator.java:139-141).
 * PercolatorOutputToPepXML: java -cp <FragPipe>/lib/* (CmdPercolator.java:270-279) with
 * <stripped>.pin, basename, the TSVs, interact-<basename>, DDA/DIA and min-prob; the script adds
 * the absolute spectra path (CmdPercolator.java:280-300).
 * Both run in the pepXML's directory (CmdPercolator.java:214, :228), FragPipe's group dir: the
 * converter writes <cwd>/<basename>.pepXML as inputfile (PercolatorOutputToPepXML.java:406-450).
 * With the task dir as work dir that is <task>/<group>, the task dir itself for the blank group.
 *
 * @return Map: percolator and converter argv, group dir, stripped basename, interact file names,
 *         FragPipe lib dir and FragPipe version
 */
def percolatorCommands(meta, pin, edited_pin, pepxml_files, percolator_cli, args, min_prob, has_crystalc, decoy_tag, write_weights, init_weights, task) {
    if (!meta.data_type) {
        error("PERCOLATOR ${meta.id}: meta.data_type is missing (FragPipe data type of the searched file)")
    }
    if (meta.group == null) {
        error("PERCOLATOR ${meta.id}: meta.group is missing (FragPipe group of the searched file)")
    }
    if (min_prob == null || min_prob.toString().isEmpty()) {
        error("PERCOLATOR ${meta.id}: min_prob is missing (workflow percolator.min-prob)")
    }
    def pepxmls = [pepxml_files].flatten()
    def stems = pepxmls.collect { f -> f.name.replaceAll(/\.[^.]+$/, '') }
    def basenames = stems.collect { s -> s.replaceAll(/_rank\d+$/, '') }.unique()
    if (basenames.size() != 1) {
        error("PERCOLATOR ${meta.id}: pepXML files ${pepxmls.collect { f -> f.name }} do not share one basename")
    }
    def basename = basenames[0]
    def stripped = has_crystalc ? basename.replaceFirst(/_c$/, '') : basename
    if (pin.name != "${stripped}.pin") {
        error("PERCOLATOR ${meta.id}: MSFragger PIN is ${pin.name}, FragPipe expects ${stripped}.pin")
    }
    def edited = [edited_pin].flatten().findAll { f -> f }
    def weights = [init_weights].flatten().findAll { f -> f }
    if (edited && edited[0].name != "${stripped}_edited.pin") {
        error("PERCOLATOR ${meta.id}: MSBooster PIN is ${edited[0].name}, FragPipe expects ${stripped}_edited.pin")
    }
    def tools = fragpipeTools(task)
    def tools_dir = tools.tools_dir
    if (!tools_dir.endsWith('/tools')) {
        error("PERCOLATOR: task.ext.fragpipe_tools_dir must be the FragPipe tools directory, got ${tools_dir}")
    }
    def fragpipe_root = tools_dir[0..<(tools_dir.length() - '/tools'.length())]
    def opts = (percolator_cli ?: '').toString().tokenize() + args.toString().tokenize()
    def only_psms = opts.contains('--only-psms')
    def percolator = ["${tools_dir}/percolator_3_7_1/linux/percolator"] + opts + ['--num-threads', "${task.cpus}"] +
        [only_psms ? '--results-psms' : '--results-peptides', "${stripped}_percolator_target_psms.tsv",
         only_psms ? '--decoy-results-psms' : '--decoy-results-peptides', "${stripped}_percolator_decoy_psms.tsv",
         '--protein-decoy-pattern', decoy_tag] +
        // first pass of Write sub mzML: Percolator saves its weights for the second pass (CmdPercolator.java:180-183)
        (write_weights ? ['--weights', "${stripped}_percolator.weights"] : []) +
        // second pass: the first pass's weights found for the _sub.mzML input (CmdPercolator.java:185-201),
        // staged in <task dir>/weights
        (weights ? ['--init-weights', "{TASK_DIR}/weights/${weights[0].name.tokenize('/').last()}", '--static', '--override'] : []) +
        [edited ? "${stripped}_edited.pin" : "${stripped}.pin"]
    def converter = ['java', '-cp', "${fragpipe_root}/lib/*", 'org.nesvilab.fragpipe.tools.percolator.PercolatorOutputToPepXML',
        "${stripped}.pin", basename, "${stripped}_percolator_target_psms.tsv", "${stripped}_percolator_decoy_psms.tsv",
        "interact-${basename}", meta.data_type == 'DDA' ? 'DDA' : 'DIA', min_prob]
    return [
        percolator      : percolator.collect { t -> t.toString() },
        converter       : converter.collect { t -> t.toString() },
        group_dir       : (meta.group == '' ? '.' : meta.group).toString(),
        stripped        : stripped,
        interact        : stems.collect { s -> "interact-${s}.pep.xml".toString() },
        fragpipe_lib    : "${fragpipe_root}/lib".toString(),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * Abacus argv as CmdPhilosopherAbacus.configure builds it (CmdPhilosopherAbacus.java:106-135):
 * `philosopher abacus <options> [--labels --plex <n>] <group dir...>`.
 *
 * @return Map: argv, group names and FragPipe version
 */
def abacusCommand(group_workspaces, group_names, protxml, abacus_cli, args, plex, task) {
    def workspaces = [group_workspaces].flatten()
    def groups = [group_names].flatten().collect { g -> g.toString() }
    if (groups.size() < 2) {
        error("PHILOSOPHERABACUS: Multi-experiment report requires more than one experiment/group (CmdPhilosopherAbacus.java:76-87), got ${groups}")
    }
    if (groups.contains('')) {
        error("PHILOSOPHERABACUS: the blank group has no directory of its own in FragPipe's work dir; abacus cannot report it")
    }
    if (workspaces.size() != groups.size()) {
        error("PHILOSOPHERABACUS: ${workspaces.size()} workspaces for groups ${groups}")
    }
    def protxmls = [protxml].flatten().findAll { f -> f }
    if (protxmls && protxmls[0].toString().tokenize('/').last() != 'combined.prot.xml') {
        error("PHILOSOPHERABACUS: Abacus reads combined.prot.xml from its directory (philosopher lib/aba/pro.go)")
    }
    def tools = fragpipeTools(task)
    def argv = [tools.philosopher, 'abacus'] + (abacus_cli ?: '').toString().tokenize() + args.toString().tokenize()
    if (plex && plex.toString().toInteger() > 0) {
        argv += ['--labels', '--plex', plex.toString()]
    }
    argv += groups
    return [
        argv            : argv.collect { t -> t.toString() },
        groups          : groups,
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * Database annotation, filter and report commands of one FragPipe group.
 *
 * {WORKDIR} marks paths under the task dir, which stands for FragPipe's work dir.
 * DbAnnotate: first group only, `database --annotate <fasta> --prefix <tag>` (CmdPhilosopherDbAnnotate.java:58-74).
 * Filter: `filter <options> --tag <tag> --pepxml <group dir> [--dbbin <first group dir>]
 * [--protxml <protxml> [--probin <first group dir>] --razor]` (CmdPhilosopherFilter.java:79-115).
 * Report: `report <flags>` (CmdPhilosopherReport.java:47-62).
 *
 * @return Map: group dir, first group dir, whether this is the first group, argv lists,
 *         pepXML names, philosopher binary and FragPipe version
 */
def philosopherFilterCommands(meta, pepxml_files, protxml, first_group_meta, first_group, filter_cli, args, dont_use_protxml, report_cli, fasta, decoy_tag, task) {
    if (meta.group == null) {
        error("PHILOSOPHERFILTER ${meta.id}: meta.group is missing (FragPipe group)")
    }
    def is_first = meta.group == first_group
    def first_meta = [first_group_meta].flatten().findAll { f -> f }
    if (is_first == (first_meta.size() > 0)) {
        error("PHILOSOPHERFILTER ${meta.id}: the first group (${first_group}) runs without first_group_meta, later groups need it")
    }
    if (meta.group == '' && !is_first) {
        error("PHILOSOPHERFILTER: the blank group sorts first in FragPipe (InputLcmsFile.compareTo)")
    }
    def protxmls = [protxml].flatten().findAll { f -> f }
    if (dont_use_protxml == (protxmls.size() > 0)) {
        error("PHILOSOPHERFILTER ${meta.id}: protXML ${protxmls} given with dont_use_protxml=${dont_use_protxml}")
    }
    def group_dir = meta.group == '' ? '.' : meta.group
    def group_path = meta.group == '' ? '{WORKDIR}' : "{WORKDIR}/${meta.group}"
    def first_path = first_group == '' ? '{WORKDIR}' : "{WORKDIR}/${first_group}"
    def tools = fragpipeTools(task)
    def philosopher = tools.philosopher
    def filter = [philosopher, 'filter'] + ((filter_cli ?: '').toString().tokenize() + args.toString().tokenize()) +
        ['--tag', decoy_tag, '--pepxml', group_path]
    if (!is_first) {
        filter += ['--dbbin', first_path]
    }
    if (!dont_use_protxml) {
        filter += ['--protxml', "{WORKDIR}/${protxmls[0].name}"]
        if (!is_first) {
            filter += ['--probin', first_path]
        }
        filter += ['--razor']
    }
    return [
        group_dir       : group_dir.toString(),
        first_dir       : first_group == '' ? '.' : first_group.toString(),
        is_first        : is_first,
        annotate        : is_first ? [philosopher, 'database', '--annotate', "{WORKDIR}/${fasta.name}", '--prefix', decoy_tag].collect { t -> t.toString() } : [],
        filter          : filter.collect { t -> t.toString() },
        report          : ([philosopher, 'report'] + (report_cli ?: '').toString().tokenize()).collect { t -> t.toString() },
        pepxml_names    : [pepxml_files].flatten().collect { f -> f.toString().tokenize('/').last() },
        philosopher     : philosopher,
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * A bash word for an argv token; {WORKDIR} (the task dir standing for FragPipe's work dir) becomes "$PWD".
 *
 * @param token Argv token
 * @return Bash word
 */
def workdirWord(token) {
    return token.startsWith('{WORKDIR}') ? "\"\$PWD\"" + shellQuote(token.substring('{WORKDIR}'.length())) : shellQuote(token)
}

/**
 * GOMEMLIMIT in GiB for Philosopher (FragPipe sets it to its RAM setting, CmdPhilosopherFilter.java:117):
 * the task memory minus a fixed JVM-free headroom for the Go runtime and the shell.
 *
 * @param task Task
 * @return GiB, at least 1
 */
def goMemLimitGb(task) {
    def headroom_gb = 2
    if (!task.memory) {
        error("PHILOSOPHERFILTER: task.memory is not set; Philosopher needs GOMEMLIMIT from the task memory")
    }
    return Math.max(1, task.memory.toGiga() - headroom_gb)
}

/**
 * Staged location of each pepXML of a multi-group step and its FragPipe location. The pepXMLs are staged
 * one per numbered directory (stageAs 'pepxml_inputs/?/*', in input order) because groups can share a
 * name: combined PeptideProphet writes interact.pep.xml in every group (CmdPeptideProphet.java:302-336).
 *
 * @param paths FragPipe locations ("<group>/<file>" or "<file>"), in the order of the staged files
 * @return [staged path, FragPipe location] per pepXML
 */
def pepxmlStagedLinks(paths) {
    return paths.withIndex().collect { rel, i -> ["pepxml_inputs/${i + 1}/${rel.tokenize('/').last()}".toString(), rel] }
}

/**
 * iProphet argv as CmdIprophet.configure builds it (CmdIprophet.java:41-58):
 * `philosopher iprophet --decoy <tag> --nonsp --output combined --threads <n> <pepxml...>`.
 *
 * @return Map: argv without the pepXML paths, the pepXML locations and the FragPipe version
 */
def iprophetCommand(pepxml_files, pepxml_paths, iprophet_cli, args, task) {
    def files = [pepxml_files].flatten().collect { f -> f.toString().tokenize('/').last() }
    def paths = [pepxml_paths].flatten().collect { p -> p.toString() }
    if (paths.collect { p -> p.tokenize('/').last() } != files || paths.unique(false).size() != paths.size()) {
        error("IPROPHET: pepXML locations ${paths} do not match the staged pepXMLs ${files}")
    }
    def tools = fragpipeTools(task)
    def argv = [tools.philosopher, 'iprophet'] + (iprophet_cli ?: '').toString().tokenize() +
        ['--threads', "${task.cpus}"] + args.toString().tokenize()
    return [
        argv            : argv.collect { t -> t.toString() },
        paths           : paths,
        links           : pepxmlStagedLinks(paths),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * PeptideProphet and RewritePepxml commands as CmdPeptideProphet builds them.
 *
 * The pepXMLs sit in FragPipe's group dir: <task>/<group>, the task dir itself for the blank group.
 * Per-file mode (CmdPeptideProphet.java:226-288): run in <group dir>/fragpipe-<pepxml>-temp with
 * `peptideprophet <opts> --decoy <tag> --database <fasta> ../<pepxml>`, output interact-<stem>.pep.xml.
 * Combine mode (CmdPeptideProphet.java:302-336): run in the group dir with
 * `peptideprophet <opts> --decoy <tag> --database <fasta> --combine <pepxml...>`, output interact.pep.xml.
 * PeptideProphet writes the interact file's absolute path (stylesheet href, summary_xml, inputfile) into it.
 * RewritePepxml: java -cp <FragPipe>/lib/* org.nesvilab.fragpipe.util.RewritePepxml <interact> <lcms...>
 * (CmdPeptideProphet.java:426-451).
 *
 * @return Map: philosopher binary, group dir, run dir, peptideprophet argv, rewrite argv prefix,
 *         interact name and its task-relative path, lcms names, run dir below the group dir and run names
 *         (peptideProphetDirLines),
 *         FragPipe lib dir and version
 */
def peptideProphetCommands(meta, pepxml_files, lcms_names, peptideprophet_cli, args, combine, fasta, decoy_tag, task) {
    if (meta.group == null) {
        error("PEPTIDEPROPHET ${meta.id}: meta.group is missing (FragPipe group of the pepXMLs)")
    }
    def group_dir = meta.group == '' ? '.' : meta.group.toString()
    def pepxmls = [pepxml_files].flatten()
    if (!combine && pepxmls.size() != 1) {
        error("PEPTIDEPROPHET ${meta.id}: per-file mode takes one pepXML, got ${pepxmls.collect { f -> f.name }}")
    }
    def lcms = [lcms_names].flatten().collect { n -> n.toString() }
    if (!lcms) {
        error("PEPTIDEPROPHET ${meta.id}: no LC-MS file names for RewritePepxml")
    }
    def tools = fragpipeTools(task)
    def tools_dir = tools.tools_dir
    if (!tools_dir.endsWith('/tools')) {
        error("PEPTIDEPROPHET: task.ext.fragpipe_tools_dir must be the FragPipe tools directory, got ${tools_dir}")
    }
    def fragpipe_root = tools_dir[0..<(tools_dir.length() - '/tools'.length())]
    // FragPipe 24.0's Philosopher binary (tools/Philosopher/philosopher-v5.1.3-RC9)
    def philosopher = tools.philosopher
    def opts = (peptideprophet_cli ?: '').toString().tokenize() + args.toString().tokenize()
    def common = [philosopher, 'peptideprophet'] + opts + ['--decoy', decoy_tag, '--database', "{CWD}/${fasta.name}"]
    def peptideprophet
    def run_dir
    def interact
    if (combine) {
        peptideprophet = common + ['--combine'] + pepxmls.collect { f -> f.name }
        run_dir = group_dir
        interact = 'interact.pep.xml'
    } else {
        def name = pepxmls[0].name
        def matched = ['pepxml', 'pep.xml'].find { ext -> name.toLowerCase().endsWith(ext) }
        if (matched == null) {
            error("PEPTIDEPROPHET ${meta.id}: could not identify the pepXML extension of ${name} (CmdPeptideProphet.java:86-99)")
        }
        peptideprophet = common + ["../${name}"]
        run_dir = group_dir == '.' ? "fragpipe-${name}-temp" : "${group_dir}/fragpipe-${name}-temp"
        interact = "interact-${name[0..<name.toLowerCase().lastIndexOf(matched)]}pep.xml"
    }
    return [
        philosopher     : philosopher,
        group_dir       : group_dir,
        run_dir         : run_dir.toString(),
        run_sub         : (combine ? '.' : "fragpipe-${pepxmls[0].name}-temp").toString(),
        peptideprophet  : peptideprophet.collect { t -> t.toString() },
        rewrite         : ['java', '-cp', "${fragpipe_root}/lib/*", 'org.nesvilab.fragpipe.util.RewritePepxml'].collect { t -> t.toString() },
        interact        : interact.toString(),
        interact_path   : (group_dir == '.' ? interact : "${group_dir}/${interact}").toString(),
        lcms            : lcms,
        run_names       : (pepxmls.collect { f -> f.name } + lcms.collect { n -> n.tokenize('/').last() }).unique(),
        fragpipe_lib    : "${fragpipe_root}/lib".toString(),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * Shell lines that set PP_GROUP_DIR, the directory PeptideProphet reads the group's pepXMLs from.
 *
 * Philosopher passes InteractParser `-a<pepXML dir>` (philosopher lib/ext/peptideprophet/peptideprophet.go:
 * 128-131, 190-193; the dir of filepath.Abs(pepXML), i.e. of the shell's logical $PWD), and InteractParser
 * joins it with the basename of each run's base_name in a 256-byte buffer: a joined path of 256
 * characters or more crashes it (SIGSEGV, Philosopher "Cannot execute program"). Measured with FragPipe
 * 24.0's Philosopher v5.1.3-RC9 on a TMT18 run file with a 48-character base name: a 206-character dir
 * passes, 207 crashes; the pepXML, output, -D database and --temp paths may be longer. FragPipe runs in
 * its work dir (CmdPeptideProphet.java:229, :332), so the group dir <task>/<group> is used as it is when
 * the path fits, which keeps interact.pep.xml as FragPipe writes it. Otherwise (long work dir, long group
 * name) PeptideProphet reads the group dir through a symlink in the task's own /tmp (mktemp -d), removed
 * when the script exits; interact.pep.xml then names that directory where it names the group dir (same
 * results, a different directory prefix, as DIA-NN's diannInputDirLines). ProteinProphet and iProphet
 * were measured fine with a 230-character dir.
 *
 * @param group_dir  directory holding the pepXMLs, relative to the task dir ('.' for the task dir)
 * @param run_names  pepXML and LC-MS file names; a run's base name is one of them without its extension,
 *                   so the full name bounds it from above
 * @return shell lines
 */
def peptideProphetDirLines(group_dir, run_names) {
    def longest = run_names.collect { n -> n.toString().length() }.max()
    def group = group_dir == '.' ? '"$PWD"' : "\"\$PWD\"/${shellQuote(group_dir)}"
    def max_path = 255
    return [
        "PP_GROUP_DIR=${group}",
        "if [ \$(( \${#PP_GROUP_DIR} + 1 + ${longest} )) -gt ${max_path} ]; then",
        '    PP_LINK_DIR="$(mktemp -d /tmp/pp.XXXXXX)"',
        '    trap \'rm -rf "$PP_LINK_DIR"\' EXIT',
        '    ln -s "$PP_GROUP_DIR" "$PP_LINK_DIR/g"',
        '    PP_GROUP_DIR="$PP_LINK_DIR/g"',
        'fi',
    ].join('\n')
}

/**
 * ProteinProphet argv as CmdProteinProphet.createCmdStub builds it (CmdProteinProphet.java:230-247).
 *
 * `philosopher proteinprophet <cmd-opts>`, then `--output combined` unless the options already
 * contain --output, then the filelist (added by the script). Output: <output>.prot.xml.
 *
 * @return Map: argv without the filelist, output file name, filelist entries, binary and versions
 */
def proteinProphetCommand(pepxml_files, filelist_paths, proteinprophet_cli, args, task) {
    def files = [pepxml_files].flatten().collect { f -> f.toString().tokenize('/').last() }
    def filelist = [filelist_paths].flatten().collect { p -> p.toString() }
    if (filelist.collect { p -> p.tokenize('/').last() } != files) {
        error("PROTEINPROPHET: filelist ${filelist} does not match the staged pepXMLs ${files}")
    }
    if (filelist.unique(false).size() != filelist.size()) {
        error("PROTEINPROPHET: duplicate pepXML in the filelist ${filelist}")
    }
    def tools = fragpipeTools(task)
    def philosopher = tools.philosopher
    def opts = (proteinprophet_cli ?: '').toString().tokenize() + args.toString().tokenize()
    def argv = [philosopher, 'proteinprophet'] + opts
    if (!argv.contains('--output')) {
        argv += ['--output', 'combined']
    }
    def output = "${argv[argv.indexOf('--output') + 1]}.prot.xml"
    return [
        argv            : argv.collect { t -> t.toString() },
        output          : output.toString(),
        filelist        : filelist,
        links           : pepxmlStagedLinks(filelist),
        philosopher     : philosopher.toString(),
        philosopher_dir : "${tools.tools_dir}/Philosopher".toString(),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * PTMProphet argv as CmdPtmProphet.configure builds it (CmdPtmProphet.java:107-124).
 *
 * Binary tools/PTMProphet/PTMProphetParser-6.3.2 (CmdPtmProphet.java:54-55, :74-79); options
 * split on whitespace with every MAXTHREADS= dropped, then MAXTHREADS=1, the pepXML name and
 * the name with .pep.xml replaced by .mod.pep.xml.
 *
 * @return Map: argv, output name, PTMProphet and FragPipe versions
 */
def ptmProphetCommand(pepxml_file, ptmprophet_cli, args, task) {
    def name = pepxml_file.toString().tokenize('/').last()
    if (!(name ==~ /.*\.pep\.xml$/)) {
        error("PTMPROPHET: ${name} is not a .pep.xml file (CmdPtmProphet.java:49 pattern1)")
    }
    def tools = fragpipeTools(task)
    def tools_dir = tools.tools_dir
    def ptmprophet_version = fragpipeToolVersions().ptmprophet
    def opts = ((ptmprophet_cli ?: '').toString().tokenize() + args.toString().tokenize()).findAll { token -> !token.startsWith('MAXTHREADS=') }
    def output = name.replaceFirst(/\.pep\.xml$/, '.mod.pep.xml')
    def argv = ["${tools_dir}/PTMProphet/PTMProphetParser-${ptmprophet_version}"] + opts + ['MAXTHREADS=1', name, output]
    def fragpipe_root = tools_dir.replaceFirst(/\/tools$/, '')
    return [
        argv              : argv.collect { t -> t.toString() },
        output            : output,
        input             : name,
        fragpipe_lib      : "${fragpipe_root}/lib".toString(),
        ptmprophet_version: ptmprophet_version,
        fragpipe_version  : tools.fragpipe_version
    ]
}

/**
 * SAINTexpress argv as CmdSaintExpress.configure builds it (CmdSaintExpress.java:84-90):
 * tools/SAINTexpress/SAINTexpress-<mode>, -R<n>, -L<n>, then cmd-opts as one element.
 *
 * @param saint_input inter.dat, bait.dat and prey.dat (CmdSaintExpress.reprintTsvToSaintInput)
 * @param mode spc or int
 * @param saint_cli "-R<max-replicates> -L<virtual-controls>"
 * @param cmd_opts saintexpress.cmd-opts, one argument when not empty
 * @param args task.ext.args
 * @param task Task (ext.fragpipe_tools_dir)
 * @return Map: argv and FragPipe version
 */
def saintexpressCommand(saint_input, mode, saint_cli, cmd_opts, args, task) {
    if (!(mode in ['spc', 'int'])) {
        error("SAINTEXPRESS: mode must be spc or int, got ${mode}")
    }
    def names = [saint_input].flatten().collect { f -> f.toString().tokenize('/').last() }.sort()
    if (names != ['bait.dat', 'inter.dat', 'prey.dat']) {
        error("SAINTEXPRESS: needs inter.dat, bait.dat and prey.dat (CmdSaintExpress.java:154-167), got ${names}")
    }
    def tools = fragpipeTools(task)
    def argv = ["${tools.tools_dir}/SAINTexpress/SAINTexpress-${mode}"] + (saint_cli ?: '').toString().tokenize()
    if (cmd_opts) {
        argv += [cmd_opts.toString()]
    }
    argv += args.toString().tokenize()
    return [
        argv            : argv.collect { t -> t.toString() },
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * One FragPipe group's TMT annotation file body ("channel sample" per line, the format
 * TmtiPanel.generateDummyAnnotationFile writes and Philosopher's --annot reads, TmtiPanel.java:712-720,
 * philosopher lib/uti/uti.go GetLabelNames) from the multi-plex TMT-Integrator annotation rows.
 *
 * @param annotation_lines Lines of the annotation TSV (plex, channel, sample, sample_name, ...; header first)
 * @param group FragPipe group ('' for the blank group)
 * @return File body, one "<channel> <sample_name>" line per row of the group, each ending in a newline
 */
def labelquantAnnotationText(annotation_lines, group) {
    def lines = [annotation_lines].flatten().collect { line -> line.toString() }
    if (!lines) {
        error("Need to specify TMT file annotations in TMT-Integrator configuration.")
    }
    def rows = lines.drop(1).collect { line -> line.split('\t', -1) }.findAll { fields -> fields[0] == group.toString() }
    if (!rows) {
        error("TMT annotation has no rows for group '${group}'")
    }
    return rows.collect { fields -> "${fields.size() > 1 ? fields[1] : ''} ${fields.size() > 3 ? fields[3] : ''}\n" }.join('').toString()
}

/**
 * FragPipe's TMT-Integrator work dir, config and argv (CmdTmtIntegrator.java:144-170, 254-300).
 *
 * The blank group's psm.tsv/protein.tsv sit at the top of the work dir (LcmsFileGroup.outputDir(wd) for group ""),
 * so TMT-Integrator names that plex after it; named groups are <work_dir>/<group>/. The psm.tsv arguments are
 * absolute, groups in FragPipe's TreeMap order (blank group first, then by name). The parser's config has
 * the output placeholder line `output: @TMT_REPORT_DIR@`, which becomes <work dir>/<report_dir>.
 * Paths are tokens under {TASK_DIR} (the task directory).
 *
 * @param group_dirs Staged group directories (input_results/<name>)
 * @param root_group_dir The group_dirs name holding the blank group ('' when none)
 * @param config_text tmt-integrator-conf.yml text from the workflow parser
 * @param work_dir Work-dir name
 * @param report_dir Output directory name (tmt-report, tmt-report-unmod)
 * @param task Task (ext.fragpipe_tools_dir)
 * @return Map: config (text with {TASK_DIR}), psm_files, stage (bash lines), jar, fragpipe_version
 */
def tmtIntegratorCommands(group_dirs, root_group_dir, config_text, work_dir, report_dir, task) {
    def placeholder = 'output: @TMT_REPORT_DIR@'
    def text = config_text.toString()
    if (!text.readLines().any { line -> line.trim() == placeholder }) {
        error("TMTINTEGRATOR: the config is not the parser's tmt-integrator-conf.yml (no '${placeholder}')")
    }
    def config = text.readLines().collect { line ->
        line.trim() == placeholder ? line.replace('@TMT_REPORT_DIR@', "{TASK_DIR}/${work_dir}/${report_dir}") : line
    }.join('\n') + (text.endsWith('\n') ? '\n' : '')
    def names = [group_dirs].flatten().collect { d -> d.toString().tokenize('/').last() }
    def root = (root_group_dir ?: '').toString()
    if (root && !names.contains(root)) {
        error("TMTINTEGRATOR: blank-group directory ${root} was not staged (${names})")
    }
    def named = names.findAll { n -> n != root }.sort()
    def psm_files = (root ? ["{TASK_DIR}/${work_dir}/psm.tsv"] : []) + named.collect { n -> "{TASK_DIR}/${work_dir}/${n}/psm.tsv" }
    if (!psm_files) {
        error("TMTINTEGRATOR: no psm.tsv staged for TMT-Integrator")
    }
    def checks = (root ? [root] : []) + named
    def stage = ["mkdir -p ${shellQuote("${work_dir}/${report_dir}")}"] +
        checks.collect { n -> "[ -f ${shellQuote("input_results/${n}/psm.tsv")} ] || { echo ${shellQuote("ERROR: group directory ${n} has no psm.tsv")} >&2; exit 1; }" } +
        (root ? ["for f in ${shellQuote("input_results/${root}")}/*; do ln -s \"\$(realpath \"\$f\")\" ${shellQuote(work_dir)}/\"\$(basename \"\$f\")\"; done"] : []) +
        named.collect { n -> "ln -s \"\$(realpath ${shellQuote("input_results/${n}")})\" ${shellQuote("${work_dir}/${n}")}" }
    def tools = fragpipeTools(task)
    return [
        config          : config.toString(),
        psm_files       : psm_files.collect { p -> p.toString() },
        stage           : stage.collect { l -> l.toString() },
        jar             : "${tools.tools_dir}/TMT-Integrator-${fragpipeToolVersions().tmtintegrator}.jar".toString(),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * A bash word for a {TASK_DIR} token: "$TASK_DIR" (set to the task directory by the module) plus the quoted rest.
 *
 * @param token Path token
 * @return Bash word
 */
def taskDirWord(token) {
    return token.startsWith('{TASK_DIR}') ? "\"\$TASK_DIR\"" + shellQuote(token.substring('{TASK_DIR}'.length())) : shellQuote(token)
}

/**
 * Bash that writes a per-task file body exactly, {TASK_DIR} expanding to "$TASK_DIR": an unquoted here-doc
 * (body escaped), through printf when the body has no final newline.
 *
 * @param name File path to write
 * @param text File body
 * @return Bash command
 */
def taskFileCommand(name, text) {
    def body = text.toString()
    def delimiter = 'TASK_FILE_END'
    if (body.readLines().contains(delimiter)) {
        error("${name} contains the here-doc delimiter line ${delimiter}")
    }
    def escaped = body.replace('\\', '\\\\').replace('$', '\\$').replace('`', '\\`')
        .replace('{TASK_DIR}', '${TASK_DIR}').with { t -> expandJarTokens(t) }
    if (body.endsWith('\n')) {
        return "cat > ${shellQuote(name)} <<${delimiter}\n${escaped.substring(0, escaped.length() - 1)}\n${delimiter}".toString()
    }
    return "printf '%s' \"\$(cat <<${delimiter}\n${escaped}\n${delimiter}\n)\" > ${shellQuote(name)}".toString()
}

/**
 * PTM-Shepherd shepherd.config and class path as FragPipe 24.0 builds them (PtmshepherdParams.createConfig,
 * CmdPtmshepherd.java:52-55, 150-165): database, one dataset line per group (in the given order), a blank line,
 * then the parser's property block with @THREADS@ and @TOOLS_DIR@ filled in. Paths are {TASK_DIR} tokens: the
 * group copies sit in groups/<dir>, the LC-MS files in spectra/.
 *
 * @param datasets [dir name, dataset name] per group (dataset01 for the blank group, PtmshepherdParams.java:86)
 * @param config_text Property block from the workflow parser
 * @param fasta_name Staged protein database name
 * @param annotation_file Staged custom annotation file (annotation/<name>), [] when annotation_file is a keyword; fills @ANNOTATION_FILE@
 * @param task Task (cpus, ext.fragpipe_tools_dir)
 * @return Map: config (text), config_arg (the config path element with its literal double quotes, CmdPtmshepherd.java:163), glyco_files (definition files to check), class_path, thermo_dir, version, fragpipe_version
 */
def ptmShepherdCommands(datasets, config_text, fasta_name, annotation_file, task) {
    def tools = fragpipeTools(task)
    def annotation = [annotation_file].flatten().findAll { f -> f }
    if (config_text.toString().contains('@ANNOTATION_FILE@') != !annotation.isEmpty()) {
        error("PTMSHEPHERD: the config's @ANNOTATION_FILE@ token and the staged annotation file disagree (${annotation})")
    }
    def properties = config_text.toString().replace('@THREADS@', "${task.cpus}").replace('@TOOLS_DIR@', tools.argv_tools_dir)
        .replace('@ANNOTATION_FILE@', annotation ? "{TASK_DIR}/${annotation[0]}" : '')
    def rows = datasets.collect { item -> [item[0].toString(), item[1].toString()] }
    def config = "database = {TASK_DIR}/${fasta_name}\n" +
        rows.collect { dir, dataset -> "dataset = ${dataset} {TASK_DIR}/groups/${dir}/psm.tsv {TASK_DIR}/spectra\n" }.join('') +
        '\n' + properties
    def lines = properties.readLines()
    // glycan definition files FragPipe takes from its tools dir (PTMSGlycanAssignPanel.java:131-134)
    def glyco_files = ['glyco_residue_list', 'glyco_mod_list', 'glyco_oxonium_list'].collect { key ->
        def line = lines.find { l -> l.startsWith("${key} = ") }
        line ? line.substring("${key} = ".length()) : ''
    }.findAll { f -> f }
    def jars = ["${tools.tools_dir}/ptmshepherd-${fragpipeToolVersions().ptmshepherd}.jar", tools.batmass_io_jar, "${tools.tools_dir}/commons-math3-3.6.1.jar",
                "${tools.tools_dir}/hipparchus-1.8/hipparchus-core-1.8.jar", "${tools.tools_dir}/hipparchus-1.8/hipparchus-stat-1.8.jar"]
    // IonQuant joins the class path only when glyco_lda_features contains "kl" (PtmshepherdParams.needsIonQuant)
    if (lines.any { l -> l ==~ /^glyco_lda_features = .*kl.*/ }) {
        jars += [tools.ionquant_jar]
    }
    return [
        config          : config.toString(),
        config_arg      : '"{TASK_DIR}/shepherd.config"',
        groups          : rows.collect { row -> row[0] },
        glyco_files     : glyco_files.collect { f -> f.toString() },
        class_path      : jars.collect { j -> j.toString() },
        thermo_dir      : tools.thermo_dir,
        glyco_mode      : lines.any { l -> l ==~ /^run_glyco_mode = true.*/ },
        version         : fragpipeToolVersions().ptmshepherd,
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * A bash word for a token holding {TASK_DIR} anywhere: the quoted pieces joined by "$TASK_DIR".
 *
 * @param token Argv token
 * @return Bash word
 */
def taskDirArg(token) {
    return token.toString().split(java.util.regex.Pattern.quote('{TASK_DIR}'), -1).collect { piece -> shellQuote(piece) }.join('"$TASK_DIR"')
}

/**
 * O-Pair's filelist_opair.txt and argv as CmdOPair builds them (CmdOPair.java:100-190): the parser's argv
 * (one token per flag, "-b 20.0") with its run-time tokens filled; FragPipe stops when a tools-dir file
 * it names is missing. Paths are {TASK_DIR} tokens; the group copy is <group>/, the spectra spectra/.
 *
 * @param group FragPipe group directory name (meta.id)
 * @param filelist LC-MS file names O-Pair reads (staged in spectra/)
 * @param argv Parser argv with @TOOLS_DIR@, @THREADS@, @FILELIST@, @PSM@, @OUTDIR@, @OXONIUM_FILE@
 * @param oxonium_file Staged user oxonium ion definitions file (oxonium/<name>) for @OXONIUM_FILE@, [] when the argv has none
 * @param task Task (cpus, ext.fragpipe_tools_dir)
 * @return Map: filelist (text), spectra names, argv, tool_files to check, dll, fragpipe_version
 */
def opairCommands(group, filelist, argv, oxonium_file, task) {
    def tools = fragpipeTools(task)
    def oxonium = [oxonium_file].flatten().findAll { f -> f }
    if ([argv].flatten().any { token -> token.toString().contains('@OXONIUM_FILE@') } != !oxonium.isEmpty()) {
        error("OPAIR: the argv's @OXONIUM_FILE@ token and the staged oxonium file disagree (${oxonium})")
    }
    def names = [filelist].flatten().collect { n -> n.toString() }
    def group_dir = "{TASK_DIR}/${group}"
    def argv_tools_dir = tools.argv_tools_dir
    def filled = [argv].flatten().collect { token ->
        token.toString().replace('@TOOLS_DIR@', argv_tools_dir).replace('@THREADS@', "${task.cpus}")
            .replace('@FILELIST@', "${group_dir}/filelist_opair.txt").replace('@PSM@', "${group_dir}/psm.tsv").replace('@OUTDIR@', group_dir)
            .replace('@OXONIUM_FILE@', oxonium ? "{TASK_DIR}/${oxonium[0]}" : '')
    }
    // glycan residue/mod definitions and the default oxonium filter come from FragPipe's tools dir (CmdOPair.java:150-190)
    def tool_files = filled.collect { token -> token.length() > 3 ? token.substring(3) : '' }.findAll { value -> value.startsWith("${argv_tools_dir}/") }
    return [
        filelist        : names.collect { n -> "{TASK_DIR}/spectra/${n}\n" }.join('').toString(),
        spectra         : names,
        argv            : filled.collect { t -> t.toString() },
        tool_files      : tool_files.collect { f -> f.toString() },
        dll             : "${tools.tools_dir}/opair/CMD.dll".toString(),
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * PairScans' argv after the JVM options as CmdPairScans builds it (CmdPairScans.java:108-127):
 * -cp <FragPipe jar>:<batmass-io> org.nesvilab.fragpipe.util.PairScans <LC-MS file> <threads>
 * <activation1> <activation2> <reverse_scan_order> <single_scan_type>. The LC-MS file is a {TASK_DIR} path;
 * PairScans writes <basename>.pairs next to it.
 *
 * @param spectra_path Staged LC-MS file, relative to the task dir (spectra/<name>)
 * @param pairscans_args The parser's four PairScans arguments
 * @param task Task (cpus, ext.fragpipe_tools_dir)
 * @return Map: argv, jars (class path entries to check), thermo_dir, bruker_dir, fragpipe_version
 */
def pairScansCommands(spectra_path, pairscans_args, task) {
    def tools = fragpipeTools(task)
    def jars = [tools.fragpipe_jar, tools.batmass_io_jar]
    def argv = ['-cp', jars.join(':'), 'org.nesvilab.fragpipe.util.PairScans', "{TASK_DIR}/${spectra_path}", "${task.cpus}"] +
        [pairscans_args].flatten()
    return [
        argv            : argv.collect { t -> t.toString() },
        jars            : jars.collect { j -> j.toString() },
        thermo_dir      : tools.thermo_dir,
        bruker_dir      : tools.bruker_dir,
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * SpecLibGen's gen_con_spec_lib.py argv, filelist and environment as CmdSpecLibGen builds them
 * (CmdSpecLibGen.java:114-216), run in the group work dir (the task dir, {TASK_DIR}).
 *
 * argv: python3 -u <tools>/speclib/gen_con_spec_lib.py <fasta> <wd> unused <wd> True unused use_easypqp
 * <rt>:<im> <threads> <convert args> <library args> <intermediate> <filelist>. A user RT/IM calibration file
 * replaces its keyword (CmdSpecLibGen.java:147-152). FRAGPIPE_PARSE_MZML_CMD is List.toString() of the java
 * command followed by the MSFragger jar path (CmdSpecLibGen.java:204-212).
 *
 * @param argv_parts Parser values: rt_cal, im_cal, convert_args, library_args, intermediate
 * @param fasta_name Staged protein database name
 * @param spectra Staged LC-MS paths (<name> in the task dir), in LC-MS table order
 * @param rt_cal_file Staged RT calibration path, or '' / [] for the keyword
 * @param im_cal_file Staged IM calibration path, or '' / [] for the keyword
 * @param ram_gb JVM heap of MSFragger's mzML parser (GB)
 * @param task Task (cpus, ext.fragpipe_tools_dir)
 * @return Map: argv, filelist (text), script, parse_mzml_cmd, fragpipe_jar, fragpipe_version
 */
def speclibgenCommand(argv_parts, fasta_name, spectra, rt_cal_file, im_cal_file, ram_gb, task) {
    def tools = fragpipeTools(task)
    def parts = argv_parts ?: [:]
    ['rt_cal', 'im_cal', 'convert_args', 'library_args', 'intermediate'].each { key ->
        if (!parts.containsKey(key)) {
            error("SPECLIBGEN: tool_configs.json speclibgen.argv_parts has no ${key} (PARSEFRAGPIPEWORKFLOW)")
        }
    }
    def rt_file = [rt_cal_file].flatten().findAll { f -> f }
    def im_file = [im_cal_file].flatten().findAll { f -> f }
    def rt = rt_file ? "{TASK_DIR}/${rt_file[0]}" : parts.rt_cal.toString()
    def im = im_file ? "{TASK_DIR}/${im_file[0]}" : parts.im_cal.toString()
    if (!rt || !im) {
        error("SPECLIBGEN: RT/IM calibration option missing (RT='${rt}', IM='${im}'): stage the user calibration file")
    }
    def files = [spectra].flatten().collect { f -> f.toString() }
    if (!files) {
        error("SPECLIBGEN: no spectra staged for SpecLibGen")
    }
    def script = "${tools.tools_dir}/speclib/gen_con_spec_lib.py".toString()
    def argv = ['python3', '-u', script, "{TASK_DIR}/${fasta_name}", '{TASK_DIR}', 'unused', '{TASK_DIR}', 'True', 'unused', 'use_easypqp',
        "${rt}:${im}", "${task.cpus}", parts.convert_args, parts.library_args, parts.intermediate, '{TASK_DIR}/filelist_speclibgen.txt']
    return [
        argv            : argv.collect { t -> t.toString() },
        filelist        : files.collect { f -> "{TASK_DIR}/${f}\n" }.join('').toString(),
        script          : script,
        parse_mzml_cmd  : "[java, -jar, -Dfile.encoding=UTF-8, -Xmx${ram_gb}G]${tools.msfragger_jar}".toString(),
        msfragger_jar   : tools.msfragger_jar,
        fragpipe_jar    : tools.fragpipe_jar,
        fragpipe_version: tools.fragpipe_version
    ]
}

/**
 * Versions of the FragPipe 24.0 tools the part-B modules run, one source for the script and the stub
 * (versions topic). Each binary is a fixed file of the FragPipe 24.0 tools dir, so its version is fixed:
 * Philosopher tools/Philosopher/philosopher-v5.1.3-RC9 reports version=v5.1.3; Percolator
 * tools/percolator_3_7_1 reports "Percolator version 3.07.1"; SAINTexpress reports "SAINT express version 3.6.1";
 * tools/TMT-Integrator-6.1.3.jar; tools/ptmshepherd-3.0.11.jar; tools/PTMProphet/PTMProphetParser-6.3.2;
 * fragpipe-speclib is pinned by FragPipe 24.0 python/python_packages/requirements.txt (fragpipe-speclib==0.1.58).
 *
 * @return Map tool -> version
 */
def fragpipeToolVersions() {
    return [
        philosopher     : '5.1.3',
        percolator      : '3.07.1',
        saintexpress    : '3.6.1',
        tmtintegrator   : '6.1.3',
        ptmshepherd     : '3.0.11',
        ptmprophet      : '6.3.2',
        fragpipe_speclib: '0.1.58',
    ]
}

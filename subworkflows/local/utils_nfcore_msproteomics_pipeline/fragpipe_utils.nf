/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FragPipe Utility Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared helper functions for FragPipe Nextflow subworkflows.
    Import via: include { shouldRunTool; getToolArgs } from './fragpipe_utils'

    These functions accept a JSON text string (from tool_configs value channel)
    and parse it on demand. JSON parsing of the ~2-5KB config is microseconds,
    so the overhead of repeated parsing is negligible.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

/**
 * Check if a tool should run based on the tool_configs JSON string.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool to check (e.g., 'msfragger', 'msbooster')
 * @return true if the tool should run, false otherwise
 */
def shouldRunTool(json_text, tool_name) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.run ?: false
}

/**
 * Get the args string for a tool from the tool_configs JSON string.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Args string for the tool, or empty string if not found
 */
def getToolArgs(json_text, tool_name) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.args ?: ''
}

/**
 * Get the IonQuant mod masses (modmasses_ionquant.txt lines) of a tool entry.
 *
 * FragPipe writes modmasses_ionquant.txt and passes --modlist for every IonQuant run, even
 * with no masses, except when O-Pair runs (FragpipeRun.java:1885; the parser then emits null).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool ('ionquant', 'ionquant_ms1', 'ionquant_isobaric')
 * @return List of mass strings (possibly empty), or false when IonQuant gets no --modlist
 */
def getToolModmasses(json_text, tool_name) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    def masses = configs?.get(tool_name)?.modmasses
    return masses instanceof List ? masses : false
}

/**
 * Get an arbitrary field value for a tool from the tool_configs JSON string.
 *
 * Used to extract individual algorithm parameters that the parser emits as
 * separate JSON fields (e.g., FPOP's region_size, control_label).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool (e.g., 'fpop')
 * @param field_name Name of the field to retrieve (e.g., 'region_size')
 * @param default_val Default value if the field is not found
 * @return Field value as a string, or default_val if not found
 */
def getToolField(json_text, tool_name, field_name, default_val = '') {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.get(field_name) ?: default_val
}

/**
 * DIA-Umpire Q tiers FragPipe adds to the LC-MS list (UmpirePanel.generateQ1/2/3, CmdUmpireSe.java:170-184).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @return List of enabled tiers, e.g. ['Q1', 'Q2', 'Q3']
 */
def umpireQTiers(json_text) {
    def tiers = new groovy.json.JsonSlurper().parseText(json_text)?.diaumpire?.q_tiers
    if (tiers == null) {
        error("tool_configs.json has no diaumpire.q_tiers - PARSEFRAGPIPEWORKFLOW must emit DIA-Umpire's Q1/Q2/Q3 settings")
    }
    return tiers
}

/**
 * Get the workflow decoy tag (`database.decoy-tag`) from the tool_configs JSON string.
 *
 * Single source of the decoy prefix for every FragPipe tool, replicating
 * FragpipeRun.java:1442 (one decoyTag passed to all Cmd*.configure calls).
 * The parser resolves FragPipe's own default and rejects an empty tag, so a
 * missing field here means a malformed tool_configs.json and is a hard error.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @return Decoy tag string (e.g., 'rev_')
 */
def getDecoyTag(json_text) {
    def decoy_tag = getToolField(json_text, 'database', 'decoy_tag')
    // an empty tag stops FragPipe (FragpipeRun.java:1444-1447; headless logs it without the HTML, SwingUtils.java:1110-1112)
    if (!decoy_tag.trim()) {
        error("Decoy protein prefix is empty. Please fix it in Database tab.")
    }
    return decoy_tag
}

/**
 * Get the report_args string for a tool from the tool_configs JSON string.
 *
 * Used by philosopher filter to pass report flags (--decoys, --removecontam)
 * to the `philosopher report` command.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Report args string, or empty string if not found
 */
def getToolReportArgs(json_text, tool_name) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.report_args ?: ''
}

/**
 * Build a FragPipe headless manifest (--manifest) from unique spectra files.
 * One line per file: path, experiment, bioreplicate, data type
 * (TabWorkflow.java:1255-1282 reads the same four tab-separated columns).
 *
 * @param entries List of [filename, experiment, bioreplicate, data_type] (experiment/bioreplicate may be null)
 * @return Tab-separated manifest string
 */
def generateFragpipeManifest(entries) {
    return entries.collect { entry ->
        def experiment = entry[1] == null ? '' : entry[1]
        def bioreplicate = entry[2] == null ? '' : entry[2]
        "${entry[0]}\t${experiment}\t${bioreplicate}\t${entry[3]}"
    }.join('\n')
}

/**
 * FragPipe's LC-MS data types, as the manifest's 4th column names them (InputLcmsFile.java:53-76).
 *
 * @return The data type strings
 */
def fragpipeDataTypes() {
    return ['DDA', 'DDA+', 'DIA', 'DIA-Quant', 'DIA-Lib', 'GPF-DIA']
}

/**
 * Data type of one LC-MS file (InputLcmsFile constructor, InputLcmsFile.java:53-76).
 * A samplesheet value is one of fragpipeDataTypes() (assets/schema_input.json enforces it);
 * a blank value is guessed from the path as FragPipe does for a manifest line without a type.
 *
 * @param data_type Samplesheet data_type (null when blank)
 * @param spectra LC-MS file path as given in the samplesheet
 * @return FragPipe data type
 */
def fragpipeDataType(data_type, spectra) {
    if (data_type != null) {
        if (!fragpipeDataTypes().contains(data_type)) {
            error("Unknown data type '${data_type}' for ${spectra}; FragPipe data types: ${fragpipeDataTypes().join(', ')}")
        }
        return data_type
    }
    return guessFragpipeDataType(spectra)
}

/**
 * InputLcmsFile.guessDataType (InputLcmsFile.java:86-101), on the absolute normalised path.
 *
 * @param spectra LC-MS file path
 * @return 'DDA', 'DIA' or 'DDA+'
 */
def guessFragpipeDataType(spectra) {
    def name = spectra.toAbsolutePath().normalize().toString()
    if (name.toLowerCase().contains('dda') || ['_Q1.', '_Q2.', '_Q3.', '_diatracer.'].any { tag -> name.contains(tag) }) {
        return 'DDA'
    }
    if (name.contains('DIA')) {
        return 'DIA'
    }
    if (name.contains('WWA')) {
        return 'DDA+'
    }
    return 'DDA'
}

/**
 * MSFragger params set a file is searched with (CmdMsfragger.java:541-552): DIA and DIA-Lib
 * share the DIA params; DIA-Quant files never reach MSFragger (FragpipeRun.java:1273-1280).
 *
 * @param data_type FragPipe data type of the file
 * @return 'DDA', 'DIA', 'GPF-DIA' or 'DDA+'
 */
def fraggerParamsType(data_type) {
    if (data_type == 'DDA' || data_type == 'GPF-DIA' || data_type == 'DDA+') {
        return data_type
    }
    if (data_type == 'DIA' || data_type == 'DIA-Lib') {
        return 'DIA'
    }
    error("Unknown data type: ${data_type}")
}

/**
 * FragPipe's second pass of Write sub mzML: a <stem>_sub.mzML input with <stem>_percolator.weights below
 * its directory. FragPipe walks the input's directory for that name, for MSFragger
 * (CmdMsfragger.java:612-626) and for Percolator (CmdPercolator.java:185-201):
 *
 *     String ss = inputLcmsFile.getPath().getFileName().toString().replaceFirst("_sub\\.mzML$", "_percolator.weights");
 *     List<Path> pList = Files.walk(inputLcmsFile.getPath().toAbsolutePath().getParent()).filter(p -> {
 *       String s = p.getFileName().toString();
 *       return s.endsWith("_percolator.weights") && s.contentEquals(ss);
 *     }).collect(Collectors.toList());
 *
 * MSFragger adds --no-ptm-score-in-pin to a command when any of its LC-MS files has a match
 * (CmdMsfragger.java:613-630; one command per data-type params set, :539-560). Percolator passes
 * --init-weights <weights> --static --override for one match and stops the run for more
 * (CmdPercolator.java:192-201).
 *
 * @param spectra_files [meta, path] per LC-MS file of the run
 * The walk runs once, at the start of the run, on the head node: it reads directory listings only, and has
 * no depth limit, as Files.walk (a bounded walk could miss a weights file FragPipe finds).
 *
 * @return the same items; meta.sub_mzml_weights lists the matches (absolute Paths, sorted) of a _sub.mzML
 *         file that has any, and meta.no_ptm_score_in_pin is true on every file of a params set with one
 */
def subMzmlWeights(spectra_files) {
    def found = spectra_files.collect { _meta, spectra ->
        def name = spectra.name
        if (!name.endsWith('_sub.mzML')) {
            return []
        }
        def wanted = name.replaceFirst(/_sub\.mzML$/, '_percolator.weights')
        def matches = []
        java.nio.file.Files.walk(spectra.toAbsolutePath().getParent()).withCloseable { stream ->
            matches = stream.collect(java.util.stream.Collectors.toList())
                .findAll { p -> p.getFileName() != null && p.getFileName().toString() == wanted }
                .collect { p -> p.toAbsolutePath().normalize() }
                .sort { p -> p.toString() }
        }
        matches
    }
    def flagged = [spectra_files, found].transpose()
        .findAll { item, matches -> matches && item[0].data_type != 'DIA-Quant' }
        .collect { item, _matches -> fraggerParamsType(item[0].data_type) }
        .unique()
    return [spectra_files, found].transpose().collect { item, matches ->
        def meta = item[0]
        def extra = [:]
        if (matches) {
            extra.sub_mzml_weights = matches
        }
        if (meta.data_type != 'DIA-Quant' && fraggerParamsType(meta.data_type) in flagged) {
            extra.no_ptm_score_in_pin = true
        }
        [extra ? meta + extra : meta, item[1]]
    }
}

/**
 * Percolator's --init-weights file of a second-pass _sub.mzML input (CmdPercolator.java:185-201):
 *
 *     if (pList.size() > 1) {
 *       throw new IllegalStateException("Found more than one weights file: " + ...joining(", "));
 *     } else if (pList.size() == 1) { cmdPp.add("--init-weights"); ... }
 *
 * FragPipe's configure catches the exception and stops the run (:202-205).
 *
 * @param meta Meta of the searched file (meta.sub_mzml_weights from subMzmlWeights)
 * @return the weights file (a Path found by subMzmlWeights), or [] when the file has none
 */
def percolatorInitWeights(meta) {
    def weights = meta.sub_mzml_weights ?: []
    if (weights.size() > 1) {
        error("Found more than one weights file: ${weights.join(', ')}")
    }
    return weights ? weights[0] : []
}

/**
 * Number of MSFragger database slices of the workflow (msfragger.misc.slice-db, written to the
 * tool_configs JSON by parse_fragpipe_workflow.py).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @return Slice count (1 = no split)
 */
def getMsfraggerNumSlices(json_text) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    def slices = configs?.get('msfragger')?.num_slices
    if (slices == null) {
        error("tool_configs.json has no msfragger.num_slices")
    }
    return slices as int
}

/**
 * calibrate_mass of one MSFragger params set, as PARSEFRAGPIPEWORKFLOW wrote it.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param params_type MSFragger params set ('DDA', 'DIA', 'GPF-DIA', 'DDA+')
 * @return calibrate_mass as an int
 */
def fraggerCalibrateMass(json_text, params_type) {
    def configs = new groovy.json.JsonSlurper()
        .setType(groovy.json.JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    def value = configs?.msfragger?.calibrate_mass?.get(params_type)
    if (value == null) {
        error("tool_configs has no msfragger calibrate_mass for the ${params_type} params set")
    }
    return value as int
}

/**
 * The one MSFragger params file MSBooster reads for the whole run (CmdMSBooster.java:134-149),
 * chosen from the manifest's data types (TabWorkflow.hasDataType over every file) and whether
 * DIA-Umpire / diaTracer run.
 *
 * @param data_types Data types of all samplesheet files
 * @param umpire_runs true when DIA-Umpire runs
 * @param diatracer_runs true when diaTracer runs
 * @return Params set ('DDA', 'DIA', 'GPF-DIA' or 'DDA+'), see fraggerParamsType
 */
def msboosterFraggerParamsType(data_types, umpire_runs, diatracer_runs) {
    if (data_types.contains('DDA') || umpire_runs || diatracer_runs) {
        return 'DDA'
    }
    if (data_types.contains('DIA') || data_types.contains('DIA-Lib')) {
        return 'DIA'
    }
    if (data_types.contains('GPF-DIA')) {
        return 'GPF-DIA'
    }
    if (data_types.contains('DDA+')) {
        return 'DDA+'
    }
    error("There are no DDA, DIA, GPF-DIA, GPF-Lib, or DDA+ data.")
}

/**
 * FragPipe's configuration checks for DIA-Umpire and diaTracer (FragpipeRun.java:1360-1366,
 * :1397-1402), over the files FragPipe processes (DIA-Quant excluded, FragpipeRun.java:1273-1280).
 *
 * @param spectra_files List of [meta, spectra] (meta.data_type set)
 * @param umpire_runs true when DIA-Umpire runs
 * @param diatracer_runs true when diaTracer runs
 * @return spectra_files unchanged
 */
def checkDiaPreprocessing(spectra_files, umpire_runs, diatracer_runs) {
    def processed = spectra_files.findAll { item -> item[0].data_type != 'DIA-Quant' }
    if (umpire_runs && processed.any { item -> item[1].toString().toLowerCase().endsWith('.d') }) {
        error("DIA-Umpire is enabled but there are timsTOF (.d) files. Please double check the parameter settings.")
    }
    if (diatracer_runs && processed.any { item -> !item[1].toString().toLowerCase().endsWith('.d') && item[0].data_type in ['DIA', 'DIA-Lib', 'GPF-DIA'] }) {
        error("diaTracer is enabled but there are DIA files not in .d format. Please double check the parameter settings.")
    }
    return spectra_files
}

/**
 * FragPipe's timsTOF compatibility check (FragpipeRun.java:2420-2443): a headless run aborts
 * when Crystal-C, FreeQuant or LabelQuant would run on timsTOF input. The input is timsTOF when
 * any LC-MS file is a .d; otherwise regular MS when any is a .raw; otherwise the workflow's
 * input.data-type.im-ms radio decides (TabWorkflow.getInputDataType, TabWorkflow.java:521-539).
 * Names are each command's getCmdName() in graph order: Crystal-C (CmdCrystalc.java:54),
 * FreeQuant (CmdFreequant.java:43; also the TMT FreeQuant that LabelQuant needs,
 * FragpipeRun.java:2470) and Quant (Isobaric) (CmdLabelquant.java:50).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param file_names Names of the LC-MS input files
 * @param im_ms true when the workflow sets workflow.input.data-type.im-ms=true
 * @return json_text unchanged
 */
def checkTimsTofTools(json_text, file_names, im_ms) {
    def names = file_names.collect { name -> name.toString().toLowerCase() }
    def timstof = names.any { name -> name.endsWith('.d') } || (!names.any { name -> name.endsWith('.raw') } && im_ms)
    if (timstof) {
        def labelquant = shouldRunTool(json_text, 'labelquant')
        def incompatible = [
            [shouldRunTool(json_text, 'crystalc'), 'Crystal-C'],
            [shouldRunTool(json_text, 'freequant') || labelquant, 'FreeQuant'],
            [labelquant, 'Quant (Isobaric)']
        ].findAll { item -> item[0] }.collect { item -> item[1] }
        if (incompatible) {
            error("timsTOF data is currently not compatible with some of tools to be run: ${incompatible.join(', ')}")
        }
    }
    return json_text
}

/**
 * Meta of a DIA-Umpire / diaTracer pseudo-MS/MS file. FragPipe re-types it DDA with the source
 * file's group as its experiment and the source replicate
 * (new InputLcmsFile(out, f.getGroup(), f.getReplicate(), "DDA"), CmdUmpireSe.java:73,
 * CmdDiaTracer.java:139), so its group is getGroup() of that new file.
 *
 * @param meta Meta of the source DIA file
 * @param id Base name of the pseudo-MS/MS file
 * @param single_group true when SpecLibGen or DIA-NN runs (see speclibOrDiannRuns)
 * @return Meta of the pseudo-MS/MS file
 */
def ddaPseudoSpectraMeta(meta, id, single_group) {
    def experiment = meta.group.trim().replaceAll(/[^A-Za-z0-9_]/, '_')
    return meta + [
        id:         id,
        experiment: experiment,
        group:      fragpipeGroup(experiment ?: null, meta.bioreplicate, single_group),
        data_type:  'DDA'
    ]
}

/**
 * FragPipe group of one LC-MS file (InputLcmsFile.java:48-51 and :131-145).
 * The experiment is trimmed and every character outside [A-Za-z0-9_] becomes '_'.
 * speclib/DIA-NN run -> ""; bioreplicate set -> experiment_bioreplicate
 * (exp_bioreplicate when experiment is blank); else experiment ("" when blank).
 *
 * @param experiment Samplesheet experiment value (null when blank)
 * @param bioreplicate Samplesheet bioreplicate value (Integer, null when blank)
 * @param speclib_or_diann true when SpecLibGen or DIA-NN runs (FragPipe then uses one group)
 * @return FragPipe group name ("" is the single default group)
 */
def fragpipeGroup(experiment, bioreplicate, speclib_or_diann) {
    if (speclib_or_diann) {
        return ''
    }
    def exp = fragpipeExperiment(experiment)
    if (bioreplicate != null) {
        return exp ? "${exp}_${bioreplicate}".toString() : "exp_${bioreplicate}".toString()
    }
    return exp
}

/**
 * Experiment of an LC-MS file as InputLcmsFile stores it (InputLcmsFile.java:48-51): trimmed,
 * every character outside [A-Za-z0-9_] replaced by '_', "" when blank.
 *
 * @param experiment Samplesheet experiment value (null when blank)
 * @return FragPipe experiment name
 */
def fragpipeExperiment(experiment) {
    return experiment == null ? '' : experiment.trim().replaceAll(/[^A-Za-z0-9_]/, '_')
}

/**
 * Whether FragPipe puts every LC-MS file in the single group "":
 * InputLcmsFile.getGroup() returns "" when SpeclibPanel.isRun() or DiannPanel.isRun()
 * (InputLcmsFile.java:131-136). The panels' run checkboxes are the workflow keys
 * speclibgen.run-speclibgen (SpeclibPanel.java:136-137) and diann.run-dia-nn.
 *
 * @param workflow_props FragPipe workflow properties (fragpipeWorkflowProperties)
 * @return true when SpecLibGen or DIA-NN runs
 */
def speclibOrDiannRuns(workflow_props) {
    return workflow_props.getProperty('speclibgen.run-speclibgen') == 'true' ||
        workflow_props.getProperty('diann.run-dia-nn') == 'true'
}

/**
 * Nextflow-internal directory / task name for a FragPipe group.
 * FragPipe writes the blank group "" into the work-dir root (LcmsFileGroup.outputDir:
 * workDir.resolve("")), which a Nextflow task cannot emit as a named directory. The
 * blank group is therefore staged as 'workdir-root'. The '-' cannot occur in a real
 * group (FragPipe replaces it with '_', InputLcmsFile.java:51), so it never collides.
 * Tools never see this name: IonQuant reads the blank group from its work-dir root and
 * publishing strips it (conf/modules.config).
 *
 * @param group FragPipe group name
 * @return Non-empty directory name for the group
 */
def fragpipeGroupDir(group) {
    if (group == null) {
        error("meta.group is missing: every spectra file needs its FragPipe group (parseSamplesheet sets it)")
    }
    return group == '' ? 'workdir-root' : group
}

/**
 * Staging directory of the blank group among per-group results, for IonQuant's root group.
 *
 * @param metas Metas of the per-group results (meta.group, meta.id)
 * @return The blank group's directory name, or '' when every group is named
 */
def fragpipeRootGroupDir(metas) {
    return metas.any { meta -> meta.group == '' } ? fragpipeGroupDir('') : ''
}

/**
 * Sub-directories of an IonQuant output directory that hold a group's psm.tsv.
 * Uses NIO Files.list() (toFile() fails on S3-backed paths).
 *
 * @param dir Directory to list
 * @return Sorted list of group directories
 */
def psmGroupDirs(dir) {
    def subdirs = []
    java.nio.file.Files.list(dir).withCloseable { stream ->
        subdirs = stream.collect(java.util.stream.Collectors.toList())
            .findAll { sub -> java.nio.file.Files.isDirectory(sub) && java.nio.file.Files.exists(sub.resolve('psm.tsv')) }
            .sort { sub -> sub.toString() }
    }
    return subdirs
}

/**
 * Results directory of a Philosopher filter/report task, from its result files.
 * The modules emit the files (not the directory) so publishDir sees each file under its
 * task-relative path and can put the blank group's files at the top (conf/modules.config);
 * the directory itself, including the .meta/ workspace, is their common parent.
 *
 * @param files One result file or a list of them (all in the same directory)
 * @return The directory holding the files
 */
def groupResultsDir(files) {
    def dirs = [files].flatten().collect { f -> f.parent }.unique()
    if (dirs.size() != 1) {
        error("Philosopher results span several directories: ${dirs}")
    }
    return dirs[0]
}

/**
 * Parse the samplesheet with nf-schema (assets/schema_input.json) into one
 * [meta, spectra] item per samplesheet row. Blank optional cells (nf-schema
 * returns []) become null. meta.group follows FragPipe's getGroup() rule.
 *
 * @param input Samplesheet path (params.input)
 * @param schema Path to assets/schema_input.json
 * @param speclib_or_diann true when SpecLibGen or DIA-NN runs (see fragpipeGroup)
 * @return List of [meta, spectra] per samplesheet row
 */
def parseSamplesheet(input, schema, speclib_or_diann) {
    return samplesheetToList(input, schema).collect { row ->
        def cells = row[0].collectEntries { key, value -> [key, value == [] ? null : value] }
        def meta = [
            id:           cells.fraction ? "${cells.sample}_f${cells.fraction}".toString() : cells.sample,
            sample:       cells.sample,
            condition:    cells.condition == null ? cells.sample : cells.condition,
            label:        cells.label,
            fraction:     cells.fraction ? cells.fraction : 1,
            experiment:   cells.experiment,
            bioreplicate: cells.bioreplicate,
            group:        fragpipeGroup(cells.experiment, cells.bioreplicate, speclib_or_diann),
            data_type:    fragpipeDataType(cells.data_type, row[1])
        ]
        [meta, row[1]]
    }
}

/**
 * One [meta, spectra] per unique spectra file (FragpipeRun.java:1258-1262 takes each
 * LC-MS file once). TMT per-channel samplesheets have one row per channel; the
 * file-level keys (group, experiment, bioreplicate, fraction, data_type) must agree across those
 * rows. meta.id is the spectra basename so converted files do not collide.
 *
 * @param rows List of [meta, spectra] from parseSamplesheet
 * @return List of [meta, spectra], one per spectra file
 */
def uniqueSpectraFiles(rows) {
    def files = rows.groupBy { row -> row[1].toString() }.collect { path, file_rows ->
        def file_metas = file_rows.collect { row -> row[0].subMap(['group', 'experiment', 'bioreplicate', 'fraction', 'data_type']) }.unique()
        if (file_metas.size() != 1) {
            error("Samplesheet rows for ${path} disagree on experiment/bioreplicate/fraction/data_type: ${file_metas}")
        }
        [[id: file_rows[0][1].baseName] + file_metas[0], file_rows[0][1]]
    }
    def duplicate_ids = files.groupBy { item -> item[0].id }.findAll { _id, items -> items.size() > 1 }.keySet()
    if (duplicate_ids) {
        error("Different spectra files share a basename (FragPipe output names would collide): ${duplicate_ids}")
    }
    return files
}

/**
 * Read a FragPipe .workflow file (Java properties, as FragPipe loads it).
 *
 * @param workflow_file Path to the .workflow file
 * @return java.util.Properties of the workflow
 */
def fragpipeWorkflowProperties(workflow_file) {
    def props = new java.util.Properties()
    workflow_file.withReader { reader -> props.load(reader) }
    return props
}

/**
 * Split a QuantLabel.java reagent string ("126, 127N, ...") the way FragPipe does ([,\s]+).
 *
 * @param names Comma/space separated reagent names
 * @return List of reagent names
 */
def reagentList(names) {
    return names.split(/[,\s]+/) as List
}

/**
 * Reagent (channel) names per FragPipe label type, verbatim from QuantLabel.java:51-69.
 *
 * @return Map of label type name to its ordered reagent names
 */
def fragpipeQuantLabels() {
    return [
        'TMT-0'     : ['126'],
        'TMT-2'     : ['126', '127C'],
        'TMT-6'     : reagentList('126, 127N, 128C, 129N, 130C, 131N'),
        'TMT-10'    : reagentList('126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N'),
        'TMT-11'    : reagentList('126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C'),
        'TMT-16'    : reagentList('126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N'),
        'TMT-18'    : reagentList('126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N'),
        'TMT-35'    : reagentList('126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N, 127D, 128ND, 128CD, 129ND, 129CD, 130ND, 130CD, 131ND, 131CD, 132ND, 132CD, 133ND, 133CD, 134ND, 134CD, 135ND, 135CD'),
        'iodoTMT-6' : reagentList('126, 127, 128, 129, 130, 131'),
        'iTRAQ-4'   : reagentList('114, 115, 116, 117'),
        'iTRAQ-8'   : reagentList('113, 114, 115, 116, 117, 118, 119, 121'),
        'sCLIP-6'   : reagentList('300, 301N, 301C, 302N, 302O, 302C'),
        'IBT-16'    : reagentList('114, 115N, 115C, 116N, 116C, 117N, 117C, 118N, 118C, 119N, 119C, 120N, 120C, 121N, 121C, 122'),
        'DiLeu-12'  : reagentList('115a, 115b, 116a, 116b, 116c, 117a, 117b, 117c, 118a, 118b, 118c, 118d'),
        'DiLeu-1'   : ['114'],
        'DeAla-13'  : ['100', '101a', '101b', '102a', '102b', '102c', '103a', '103b', '103c', '104a', '104b', '105a', '105b'],
        'Trp-2'     : ['597', '603'],
        'develop'   : (1..32).collect { n -> "xtag${n}".toString() },
        'develop2'  : reagentList('114, 115a, 115b, 115c, 116a, 116b, 116c, 116d, 116e, 117a, 117b, 117c, 117d, 117e, 117f, 118a, 118b, 118c, 118d, 118e, 118f, 118g, 119a, 119b, 119c, 119d, 119e, 119f, 119g')
    ]
}

/**
 * Fail unless every spectra file is an LC-MS format FragPipe takes: .mgf, .mzML, .mzXML, .mzBIN,
 * .raw and Bruker .d (CmdMsfragger.java:764-766), none compressed. They are searched as given.
 *
 * @param spectra_files List of [meta, spectra] from uniqueSpectraFiles
 * @return spectra_files unchanged
 */
def requireFragpipeLcmsFiles(spectra_files) {
    def unsupported = spectra_files.collect { item -> item[1] }.findAll { spectra ->
        !(spectra.name.toLowerCase() ==~ /.+\.(raw|mzml|mzxml|mgf|mzbin|d)/)
    }
    if (unsupported) {
        error("FragPipe mode takes .raw, .mzML, .mzXML, .mgf, .mzBIN or Bruker .d LC-MS files (CmdMsfragger.java:764-766); compressed files are not FragPipe input. Unsupported: ${unsupported}")
    }
    return spectra_files
}

/**
 * The one LC-MS file FragPipe's CheckCentroid reads: the first of the run's files without DIA-Quant
 * (FragpipeRun.java:1273-1281), taken group by group in sorted group order (TabWorkflow.getLcmsFileGroups
 * is a TreeMap) and in manifest order within a group; CmdCheckCentroid.java:52-65 stops after one file.
 *
 * @param files List of [meta, spectra] in manifest order (meta.group, meta.data_type)
 * @return [meta, spectra] of the checked file, or null when no file is searched
 */
def checkCentroidFile(files) {
    def searched = files.findAll { item -> item[0].data_type != 'DIA-Quant' }
    if (!searched) {
        return null
    }
    def first_group = searched.collect { item -> item[0].group }.min()
    return searched.find { item -> item[0].group == first_group }
}

/**
 * MSFragger console logs in the order FragPipe runs MSFragger, the order of its console text, of
 * which the SDRF writer (TabMsfragger.java:2359-2363) and the second pass (FragpipeRun.java:598-635)
 * take the first "New ..." line. FragPipe's LC-MS list (FragpipeRun.java:1274-1281) holds the groups
 * sorted (TabWorkflow.getLcmsFileGroups is a TreeMap), manifest order within a group; DIA-Umpire and
 * diaTracer replace a file in place by its DDA-typed Q1, Q2, Q3 / _diatracer.mzML files
 * (FragpipeRun.java:1380-1382, :1421-1423; CmdUmpireSe.java:59-78, CmdDiaTracer.java:128-143).
 * CmdMsfragger.java:539-555 then runs that list bucketed in a TreeMap keyed by data type (DIA-Lib in
 * DIA; String order DDA, DDA+, DIA, GPF-DIA, not the order of the comment at :568).
 *
 * @param items List of [meta, log] (meta.id, meta.data_type of the searched file)
 * @param spectra_files List of [meta, LC-MS file] in manifest order (meta.id, meta.group)
 * @return The logs in MSFragger's run order
 */
def msfraggerConsoleLogs(items, spectra_files) {
    def buckets = ['DDA': 'DDA', 'DDA+': 'DDA+', 'DIA': 'DIA', 'DIA-Lib': 'DIA', 'GPF-DIA': 'GPF-DIA']
    // searched file id -> [group, manifest index, index of the file among the source file's outputs]
    def run_order = [:]
    spectra_files.eachWithIndex { item, index ->
        def meta = item[0]
        // ids as FRAGPIPE_CONVERT names the outputs (ddaPseudoSpectraMeta: the output file's base name)
        def name = item[1].name
        def outputs = ['_Q1', '_Q2', '_Q3'].collect { tier -> "${item[1].baseName}${tier}".toString() } +
            [name.replace('.d', '_diatracer.mzML').replaceFirst(/\.mzML$/, '')]
        run_order[meta.id] = [meta.group, index, 0]
        outputs.eachWithIndex { output_id, n -> run_order[output_id] = [meta.group, index, n + 1] }
    }
    return items
        .collect { meta, log_file ->
            if (!buckets.containsKey(meta.data_type)) {
                error("MSFragger does not search data type '${meta.data_type}' of ${meta.id} (CmdMsfragger.java:541-551)")
            }
            if (!run_order.containsKey(meta.id)) {
                error("MSFragger log of ${meta.id} is neither an LC-MS file of the manifest nor a DIA-Umpire / diaTracer output of one")
            }
            return [[buckets[meta.data_type]] + run_order[meta.id], log_file]
        }
        .sort(false) { a, b -> (a[0][0] <=> b[0][0]) ?: (a[0][1] <=> b[0][1]) ?: (a[0][2] <=> b[0][2]) ?: (a[0][3] <=> b[0][3]) }
        .collect { _key, log_file -> log_file }
}

/**
 * FragPipe's LC-MS directory after the search: each LC-MS file with the files MSFragger wrote
 * next to it (MSFRAGGER lcms_files: _uncalibrated.mzML, _calibrated.mzML, mzBIN caches).
 * Call from a workflow body; joins on meta.id.
 *
 * @param ch_spectra Channel of [meta, LC-MS file] that MSFragger searched
 * @param ch_lcms_files Channel of [meta, files] from MSFRAGGER.out.lcms_files
 * @return Channel of [meta, [LC-MS file, side files...]]
 */
def fragpipeLcmsDir(ch_spectra, ch_lcms_files) {
    return ch_spectra
        .map { meta, spectra -> [meta.id, meta, spectra] }
        .join(ch_lcms_files.map { meta, files -> [meta.id, files instanceof List ? files : [files]] }, remainder: true)
        .filter { _key, meta, _spectra, _files -> meta != null }
        .map { _key, meta, spectra, files -> [meta, [spectra] + (files ?: [])] }
}

/**
 * FragPipe's hasCalibratedMzml: write_calibrated_mzml on and mass calibration on
 * (FragpipeRun.java:1642, tabMsf.isWriteCalMzml() && tabMsf.getMassCalibration() > 0), read from
 * the MSFragger params file the workflow gives.
 *
 * @param params_text Text of the MSFragger params file (fragger.params)
 * @return true when Percolator / PeptideProphet name <base>_calibrated.mzML
 */
def fraggerWritesCalibratedMzml(params_text) {
    return fraggerParamValue(params_text, 'write_calibrated_mzml') == '1' && (fraggerParamValue(params_text, 'calibrate_mass') as int) > 0
}

/**
 * Value of a key in MSFragger params text ("key = value  # comment").
 *
 * @param params_text Text of the MSFragger params file
 * @param key Parameter name
 * @return The value, without comment; fails when the key is missing
 */
def fraggerParamValue(params_text, key) {
    def line = params_text.readLines().find { l -> l.contains('=') && l.split('=', 2)[0].trim() == key }
    if (line == null) {
        error("MSFragger params have no ${key}")
    }
    return line.split('=', 2)[1].replaceAll(/#.*/, '').trim()
}

/**
 * The spectra file Percolator and PeptideProphet name in their outputs (CmdPercolator.java:219-224,
 * CmdPeptideProphet.java:442-447): <base>_calibrated.mzML when hasCalibratedMzml, else
 * <base>_uncalibrated.mzML for .raw/.d input, else the LC-MS file itself.
 *
 * @param files [LC-MS file, files MSFragger wrote next to it...] (FRAGPIPE_SEARCH lcms)
 * @param has_calibrated_mzml fraggerWritesCalibratedMzml of the run's MSFragger params
 * @return The chosen file; fails when MSFragger did not write it
 */
def fragpipeValidationLcms(files, has_calibrated_mzml) {
    def lcms = files[0]
    def base = lcms.name.take(lcms.name.lastIndexOf('.'))
    def name = has_calibrated_mzml ? "${base}_calibrated.mzML".toString()
        : lcms.name.toLowerCase() ==~ /.+\.(raw|d)/ ? "${base}_uncalibrated.mzML".toString()
        : lcms.name
    def chosen = files.find { f -> f.name == name }
    if (!chosen) {
        error("MSFragger wrote no ${name} for ${lcms.name}; FragPipe's validation tools read it (CmdPercolator.java:219-224)")
    }
    return chosen
}

/**
 * TMT annotation, one plex per FragPipe group (TmtiPanel.getAnnotations, TmtiPanel.java:684-706),
 * written in the multi-plex TMT-Integrator format
 * (plex, channel, sample, sample_name, condition, replicate); IonQuant reads columns 2 and 4.
 * Only built when TMT-Integrator runs (FragpipeRun.java:1917-1922).
 *
 * - A group with labelled samplesheet rows is a user annotation: `sample` is the channel's
 *   sample name, `NA` excludes the channel. It must list every reagent of the label type
 *   (CmdTmtIntegrator.java:197-202), each channel once.
 * - A group without labelled rows gets FragPipe's default annotation: every reagent named
 *   "<group>_<channel>" (TmtiPanel.generateDummyAnnotationFile, TmtiPanel.java:712-720),
 *   i.e. "_<channel>" for the blank group "".
 * - Ref-tag sanity checks (CmdTmtIntegrator.java:179-197), non-NA sample names unique across all
 *   plexes (:203-208), and a reference sample when "Reference sample" is selected (:212-235),
 *   with FragPipe's messages (annotation files named <group>/<group>_annotation.txt).
 *
 * @param rows List of [meta, spectra] from parseSamplesheet
 * @param spectra_files List of [meta, spectra] from uniqueSpectraFiles (defines the groups)
 * @param workflow_props FragPipe workflow properties (fragpipeWorkflowProperties)
 * @return Annotation lines (header first), or [] when TMT-Integrator does not run
 */
def tmtAnnotationLines(rows, spectra_files, workflow_props) {
    if (workflow_props.getProperty('tmtintegrator.run-tmtintegrator') != 'true') {
        return []
    }
    // Label combo (TmtiPanel.java:323-324, 665-673): name, else the config number of
    // ui.translations:83-101 (first match), else the first item TMT-0
    def channel_value = (workflow_props.getProperty('tmtintegrator.channel_num') ?: '').trim()
    def numeric_label = ['0': 'TMT-0', '2': 'TMT-2', '6': 'TMT-6', '10': 'TMT-10', '11': 'TMT-11', '16': 'TMT-16', '18': 'TMT-18', '35': 'TMT-35',
                         '4': 'iTRAQ-4', '8': 'iTRAQ-8', '12': 'DiLeu-12', '1': 'DiLeu-1', '13': 'DeAla-13', '32': 'develop', '29': 'develop2']
    def label_type = !channel_value ? 'TMT-0' :
        (fragpipeQuantLabels().keySet().find { name -> name.equalsIgnoreCase(channel_value) } ?: numeric_label[channel_value] ?: channel_value)
    def reagents = fragpipeQuantLabels()[label_type]
    if (reagents == null) {
        error("Unknown TMT-Integrator label type tmtintegrator.channel_num=${label_type}; known: ${fragpipeQuantLabels().keySet()}")
    }
    // Sanity checks in CmdTmtIntegrator.java's order (:174-197)
    // TmtiPanel widgets (TmtiPanel.java:334-385): ref_tag starts as "Bridge", ref_d_tag as "Pool",
    // "Define reference" as "Reference sample" (-1); the checks use the config ref_tag, which
    // formToConfig cleans with QuantLabelAnnotation.unifyAnnotationSampleName
    def ref_tag = (workflow_props.getProperty('tmtintegrator.ref_tag', 'Bridge')).trim().replaceAll(/[^0-9a-zA-Z_-]/, '')
    def ref_d_tag = workflow_props.getProperty('tmtintegrator.ref_d_tag', 'Pool')
    def is_tmt35 = label_type == 'TMT-35'
    def add_ref = workflow_props.getProperty('tmtintegrator.add_Ref', '-1').trim()
    if (!(add_ref in ['-1', '1', 'Reference sample', 'Virtual'])) {
        error("No matching add_Ref value: tmtintegrator.add_Ref=${add_ref}")
    }
    def is_real_reference = add_ref in ['-1', 'Reference sample']
    if (is_tmt35 && ref_tag != ref_d_tag && (ref_tag.contains(ref_d_tag) || ref_d_tag.contains(ref_tag))) {
        error("For TMT 35-plex, the 'Ref sample tag' and 'Ref D sample tag' can't be substrings of each other.")
    }
    if (is_real_reference && !ref_tag.trim()) {
        error("'Ref sample tag' can't be empty in 'Quant (Isobaric)' tab.")
    }
    if (is_real_reference && is_tmt35 && !ref_d_tag.trim()) {
        error("'Ref D sample tag' can't be empty in 'Quant (Isobaric)' tab.")
    }
    def groups = spectra_files.collect { item -> item[0].group }.unique().sort()
    def annotations = groups.collectEntries { group ->
        def labelled = rows.findAll { row -> row[0].group == group && row[0].label }
        if (!labelled) {
            return [group, reagents.collect { channel -> [channel, "${group}_${channel}".toString(), "${group}_${channel}".toString()] }]
        }
        def by_channel = labelled.groupBy { row -> row[0].label.replaceAll(/^TMT/, '') }
        by_channel.each { channel, channel_rows ->
            if (!reagents.contains(channel)) {
                error("TMT label '${channel_rows[0][0].label}' of group '${group}' is not a channel of label type ${label_type}: ${reagents}")
            }
            def assigned = channel_rows.collect { row -> [row[0].sample, row[0].condition] }.unique()
            if (assigned.size() != 1) {
                error("Group '${group}' channel ${channel} is assigned to more than one sample/condition: ${assigned}")
            }
        }
        if (by_channel.size() != reagents.size()) {
            error("Number of the samples in the annotation file does not match the number of channels in the 'label type' of the Quant (Isobaric) tab. " +
                "Group '${group}' lists channels ${by_channel.keySet().sort()} but ${label_type} has ${reagents}; use NA as the sample name to exclude a channel.")
        }
        [group, reagents.collect { channel ->
            def meta = by_channel[channel][0][0]
            [channel, meta.sample, meta.condition]
        }]
    }
    def named = annotations.values().sum().collect { entry -> entry[1] }.findAll { sample -> !sample.equalsIgnoreCase('na') }
    def duplicated = named.countBy { sample -> sample }.findAll { _sample, n -> n > 1 }.keySet()
    if (duplicated) {
        error("Duplicate samples found in annotation files. The sample names must be unique among all experiments. Duplicates: ${duplicated}")
    }
    if (is_real_reference) {
        def missing_ref = annotations.findAll { _group, entries ->
            !entries.any { entry -> entry[1].contains(ref_tag) } ||
                (is_tmt35 && !entries.any { entry -> entry[1].contains(ref_d_tag) })
        }.keySet().collect { group -> "${group}/${group}_annotation.txt".toString() }
        if (missing_ref) {
            error("Found annotation files without reference channel\nspecified:\n" + missing_ref.join('\n') +
                "\n\nOne sample name in each annotation file must start with\nthe reference tag you set. You can change that in Quant tab, [Ref Sample Tag] text field.")
        }
    }
    def lines = annotations.collectMany { group, entries ->
        entries.collect { entry -> [group, entry[0], entry[1], entry[1], entry[2], 1].join('\t') }
    }
    return ['plex\tchannel\tsample\tsample_name\tcondition\treplicate'] + lines
}

/**
 * A tool field of tool_configs.json that must be present (false and '' are valid values).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @param field_name Field the parser emits for that tool
 * @return The field value
 */
def requireToolField(json_text, tool_name, field_name) {
    def entry = new groovy.json.JsonSlurper().parseText(json_text)?.get(tool_name)
    if (entry == null || !entry.containsKey(field_name)) {
        error("tool_configs.json has no ${tool_name}.${field_name} - PARSEFRAGPIPEWORKFLOW must emit it")
    }
    return entry[field_name]
}

/**
 * Sort key of FragPipe's LC-MS file order (InputLcmsFile.compareTo, InputLcmsFile.java:116-119):
 * getGroup2() (the group without the SpecLib/DIA-NN rule), then the data type, then the file path.
 * Spectra files sit in one directory, so the path order is the order of the file names.
 *
 * @param meta Meta with experiment, bioreplicate, data_type and spectra_name (set by FRAGPIPE_VALIDATE)
 * @return [group2, data type, spectra file name]
 */
def fragpipeLcmsOrderKey(meta) {
    if (!meta.spectra_name) {
        error("meta.spectra_name is missing for ${meta.id}: FragPipe orders LC-MS files by their path (InputLcmsFile.compareTo)")
    }
    return [fragpipeGroup(meta.experiment, meta.bioreplicate, false), meta.data_type ?: '', meta.spectra_name.toString()]
}

/**
 * Crystal-C runs numbered as FragPipe numbers their params files: one index over the whole run,
 * ++index over the TreeMap of LC-MS files in InputLcmsFile order (FragpipeRun.java:1459) and,
 * within a file, its pepXMLs by rank (CmdCrystalc.java:176-185).
 *
 * @param items [meta (with spectra_name), pepXML, LC-MS file, rank] for every Crystal-C run
 * @return [meta + crystalc_index, pepXML, LC-MS file] in FragPipe's order
 */
def crystalcRunIndex(items) {
    return items
        .sort(false) { a, b -> compareFragpipeLcmsKeys(fragpipeLcmsOrderKey(a[0]), fragpipeLcmsOrderKey(b[0])) ?: a[3] <=> b[3] }
        .withIndex()
        .collect { item, index -> [item[0] + [crystalc_index: index], item[1], item[2]] }
}

/**
 * Comparator of fragpipeLcmsOrderKey keys (String.compareTo element by element, as Java's comparator chain).
 *
 * @param a First key
 * @param b Second key
 * @return negative, zero or positive
 */
def compareFragpipeLcmsKeys(a, b) {
    def i = [0, 1, 2].find { n -> a[n] != b[n] }
    return i == null ? 0 : a[i].compareTo(b[i])
}

/**
 * Rank order of MSFragger's ranked pepXMLs (<stem>_rank<N>.<ext>, N = 1, 2, ...); other files keep their name order.
 *
 * @param files pepXML / interact files of one LC-MS file
 * @return Files sorted by rank number
 */
def sortByRank(files) {
    return [files].flatten().sort { a, b ->
        def ra = (a.name =~ /_rank(\d+)\./)
        def rb = (b.name =~ /_rank(\d+)\./)
        (ra.find() && rb.find()) ? (ra.group(1).toInteger() <=> rb.group(1).toInteger()) : (a.name <=> b.name)
    }
}

/**
 * One tool's entry of the tool_configs JSON (run flag, args and the tool's extra fields).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Map of the tool's entry
 */
def getToolEntry(json_text, tool_name) {
    def entry = new groovy.json.JsonSlurper().parseText(json_text)?.get(tool_name)
    if (entry == null) {
        error("tool_configs.json has no '${tool_name}' entry")
    }
    return entry
}

/**
 * Rows of the manifest FragPipe saves as fragpipe-files.fp-manifest (TabWorkflow.manifestSave,
 * TabWorkflow.java:1241-1250): absolute LC-MS path, experiment as InputLcmsFile stores it,
 * bioreplicate, data type; in LC-MS table order.
 * The path is the one the user gave: a local path as FragPipe writes it (toAbsolutePath().normalize()),
 * an object-store path as its URI (s3://bucket/key; S3Path.toString() drops the scheme).
 * With staged_dir the path is the file's name under that task directory instead, for a task that
 * loads the manifest (TabWorkflow.manifestLoad resolves a relative path against the task dir and
 * drops files that do not exist there).
 *
 * @param spectra_files List of [meta, spectra] (meta.experiment, meta.bioreplicate, meta.data_type)
 * @param staged_dir Directory the task stages the spectra into (stageAs: '<staged_dir>/*'), or null
 * @return List of [path, experiment, bioreplicate, data_type]
 */
def fragpipeManifestEntries(spectra_files, staged_dir = null) {
    return spectra_files.collect { item ->
        def path = staged_dir == null
            ? item[1].toAbsolutePath().normalize().toUriString()
            : "${staged_dir}/${item[1].name}".toString()
        [path, fragpipeExperiment(item[0].experiment), item[0].bioreplicate, item[0].data_type]
    }
}

/**
 * Location of a group's file relative to FragPipe's work dir: "<group>/<name>", or "<name>" for the
 * blank group, which FragPipe writes into the work dir itself (LcmsFileGroup.outputDir).
 *
 * @param group FragPipe group name
 * @param name File name
 * @return Relative path
 */
def fragpipeGroupPath(group, name) {
    return group == '' ? name.toString() : "${group}/${name}".toString()
}

/**
 * Spectra path name FragPipe writes into validated pepXMLs (CmdPercolator.java:219-224,
 * CmdPeptideProphet.java:439-450): <stem>_calibrated.mzML when MSFragger writes calibrated mzML,
 * <stem>_uncalibrated.mzML for .raw/.d, else the searched file name.
 *
 * @param name Searched LC-MS file name
 * @param has_calibrated_mzml msfragger.write_calibrated_mzml && msfragger.calibrate_mass > 0
 * @return File name
 */
def fragpipeLcmsName(name, has_calibrated_mzml) {
    def stem = name.contains('.') ? name[0..<name.lastIndexOf('.')] : name
    if (has_calibrated_mzml) {
        return "${stem}_calibrated.mzML".toString()
    }
    if (name.toLowerCase().endsWith('.raw') || name.toLowerCase().endsWith('.d')) {
        return "${stem}_uncalibrated.mzML".toString()
    }
    return name.toString()
}

/**
 * Path of an LC-MS input file as FragPipe prints it (InputLcmsFile.getPath().toAbsolutePath()).
 *
 * @param spectra Samplesheet LC-MS file
 * @return Absolute local path, or the URI for remote files
 */
def lcmsAbsolutePath(spectra) {
    return spectra.toUri().scheme == 'file' ? spectra.toAbsolutePath().toString() : spectra.toUriString()
}

/**
 * FilenameUtils.getBaseName of an LC-MS file: its name without the last extension.
 *
 * @param spectra LC-MS file
 * @return Base name
 */
def lcmsBaseName(spectra) {
    def name = spectra.name
    def dot = name.lastIndexOf('.')
    return dot < 0 ? name : name.substring(0, dot)
}

/**
 * experiment_annotation.tsv of a label-free IonQuant run, ported from
 * ToolingUtils.generateLFQExperimentAnnotation(wd, 0) (ToolingUtils.java:222-286).
 *
 * Files are FragPipe's LC-MS list grouped by getGroup2() (TreeMap, TabWorkflow.java:1651-1662),
 * in samplesheet order within a group. GPF-DIA/DIA-Lib files are skipped, and so are non-DIA
 * files when any file is DIA or DIA-Quant (filterRun). With one group named "" the sample name
 * is the file base name minus the common prefix of all base names, split on "_"; otherwise
 * sample = sample_name = getGroup2(), condition = experiment.split("_")[0], replicate =
 * bioreplicate or 1. String.split is Java's (trailing empty strings dropped).
 *
 * @param spectra_files List of [meta, spectra] of every samplesheet LC-MS file
 * @return Lines of experiment_annotation.tsv, header first
 */
def lfqExperimentAnnotationLines(spectra_files) {
    def files = spectra_files.collect { item ->
        def meta = item[0]
        [
            path:       lcmsAbsolutePath(item[1]),
            base_name:  lcmsBaseName(item[1]),
            experiment: meta.experiment == null ? '' : meta.experiment.trim().replaceAll(/[^A-Za-z0-9_]/, '_'),
            replicate:  meta.bioreplicate,
            group2:     fragpipeGroup(meta.experiment, meta.bioreplicate, false),
            data_type:  meta.data_type
        ]
    }
    def has_dia = files.any { f -> f.data_type in ['DIA', 'DIA-Quant'] }
    def groups = files.groupBy { f -> f.group2 }
    def ordered = groups.keySet().sort().collectMany { group -> groups[group] }.findAll { f ->
        !(f.data_type in ['GPF-DIA', 'DIA-Lib'] || (has_dia && !(f.data_type in ['DIA', 'DIA-Quant'])))
    }
    def lines = ['file\tsample\tsample_name\tcondition\treplicate']
    if (groups.size() == 1 && groups.keySet()[0] == '') {
        // org.apache.commons.lang3.StringUtils.getCommonPrefix over the distinct base names
        def names = files.collect { f -> f.base_name }.unique()
        def min_length = names*.length().min()
        def first_difference = (0..<min_length).find { i -> names.collect { n -> n[i] }.unique().size() > 1 }
        def prefix_length = first_difference == null ? min_length : first_difference
        ordered.each { f ->
            def sample_name = f.base_name.substring(prefix_length)
            def tt = sample_name.split('_')
            def row
            if (tt.length == 3) {
                row = tt[2].isInteger() ? [f.path, f.path, sample_name, "${tt[0]}_${tt[1]}", tt[2].toInteger()] :
                    [f.path, f.path, sample_name, "${tt[1]}_${tt[2]}", 1]
            } else if (tt.length == 2) {
                row = tt[1].isInteger() ? [f.path, f.path, sample_name, tt[0], tt[1].toInteger()] :
                    [f.path, f.path, sample_name, tt[1], 1]
            } else {
                row = [f.path, f.path, sample_name, sample_name, 1]
            }
            lines << row.join('\t')
        }
    } else {
        ordered.each { f ->
            def parts = f.experiment.split('_')
            if (parts.length == 0) {
                error("FragPipe cannot write experiment_annotation.tsv for experiment '${f.experiment}' (it has no characters besides '_')")
            }
            lines << [f.path, f.group2, f.group2, parts[0].trim(), f.replicate == null ? 1 : f.replicate].join('\t')
        }
    }
    return lines
}

/**
 * experiment_annotation.tsv of the TMT IonQuant isobaric pass and of Philosopher LabelQuant,
 * ported from ToolingUtils.writeIsobaricQuantExperimentAnnotation (ToolingUtils.java:299-334).
 *
 * Reads each group's "channel sample" annotation (the per-plex columns 2 and 4 of the
 * TMT-Integrator annotation, see tmtAnnotationLines) in group order; NA samples are skipped and
 * condition/replicate come from splitting the sample name on "_" (3 parts: second part and
 * integer third part; 2 parts: first part and integer second part; otherwise the name and 1).
 *
 * @param annotation_lines Lines of the multi-plex TMT annotation (header first)
 * @return Lines of experiment_annotation.tsv, header first
 */
def isobaricExperimentAnnotationLines(annotation_lines) {
    def lines = ['plex\tchannel\tsample\tsample_name\tcondition\treplicate']
    def rows = annotation_lines.drop(1).findAll { line -> line.trim() }.collect { line -> line.split('\t', -1) }
    def plexes = rows.collect { cols -> cols[0] }.unique().sort()
    plexes.each { plex ->
        rows.findAll { cols -> cols[0] == plex }.each { cols ->
            def channel = cols[1].trim()
            def sample = cols[3].trim()
            if (sample.equalsIgnoreCase('na')) {
                return
            }
            def parts2 = sample.split('_')
            def condition
            def replicate
            if (parts2.length == 3) {
                condition = parts2[2].isInteger() ? parts2[1].trim() : parts2[2].trim()
                replicate = parts2[2].isInteger() ? parts2[2].toInteger() : 1
            } else if (parts2.length == 2) {
                condition = parts2[1].isInteger() ? parts2[0].trim() : parts2[1].trim()
                replicate = parts2[1].isInteger() ? parts2[1].toInteger() : 1
            } else {
                condition = sample
                replicate = 1
            }
            lines << [plex, channel, sample, sample, condition, replicate].join('\t')
        }
    }
    return lines
}

/**
 * Apply FragPipe's input data type to an IonQuant option string: --ionmobility is 1 when any
 * LC-MS input is a Bruker .d, 0 when any is .raw, else the workflow's IM-MS button value the
 * parser emitted (TabWorkflow.getInputDataType, TabWorkflow.java:521-538; CmdIonquant.java:197-198).
 *
 * @param args IonQuant options from the parser
 * @param spectra_files List of [meta, spectra] of every samplesheet LC-MS file
 * @return args with the --ionmobility value FragPipe uses
 */
def ionquantInputDataTypeArgs(args, spectra_files) {
    def names = spectra_files.collect { item -> item[1].name.toLowerCase() }
    def value = names.any { n -> n.endsWith('.d') } ? '1' : (names.any { n -> n.endsWith('.raw') } ? '0' : null)
    if (value == null) {
        return args
    }
    if (!(args =~ /--ionmobility \S+/)) {
        error("IonQuant options lack --ionmobility: ${args}")
    }
    return args.replaceFirst(/--ionmobility \S+/, "--ionmobility ${value}")
}

/**
 * MSstatsTMT_annotation.csv TMT-Integrator's configure writes when Philosopher MSstats is on,
 * ported from CmdTmtIntegrator.java:105-140: per group (name order) and annotation line (NA
 * skipped), one row per LC-MS file of the group with Fraction counting the files from 1;
 * Condition is the second "_" part of the sample name when it has exactly two parts.
 *
 * @param annotation_lines Lines of the multi-plex TMT annotation (header first)
 * @param spectra_files List of [meta, spectra] of the searched LC-MS files (meta.group)
 * @return Lines of MSstatsTMT_annotation.csv, header first
 */
def msstatsTmtAnnotationLines(annotation_lines, spectra_files) {
    def lines = ['Run,Fraction,TechRepMixture,Mixture,Channel,BioReplicate,Condition']
    def rows = annotation_lines.drop(1).findAll { line -> line.trim() }.collect { line -> line.split('\t', -1) }
    rows.collect { cols -> cols[0] }.unique().sort().each { group ->
        def runs = spectra_files.findAll { item -> item[0].group == group }.collect { item -> lcmsBaseName(item[1]) }
        rows.findAll { cols -> cols[0] == group }.each { cols ->
            def channel = cols[1].trim()
            def sample = cols[3].trim()
            if (sample.equalsIgnoreCase('na')) {
                return
            }
            def parts2 = sample.split('_')
            def condition = parts2.length == 2 ? parts2[1].trim() : ''
            runs.eachWithIndex { run, i ->
                lines << [run, i + 1, 1, group, channel, sample, condition].join(',')
            }
        }
    }
    return lines
}

/**
 * Group directories an IonQuant pass leaves for the next one: the named groups inside its
 * results_dir and the blank group's copy in its root_group output.
 *
 * @param dir IonQuant results_dir
 * @param root_dir IonQuant root_group directory (ionquant_root), null when the run has no blank group
 *        (IonQuant emits it only then)
 * @return [ [group dirs], root_group_dir ] with root_group_dir '' when there is no blank group
 */
def ionquantGroupDirs(dir, root_dir) {
    def roots = root_dir == null ? [] : psmGroupDirs(root_dir)
    return [psmGroupDirs(dir) + roots, roots ? roots[0].name : '']
}

/**
 * experiment_annotation.tsv that CmdDiann writes for DIA-NN runs
 * (generateLFQExperimentAnnotation(wd, 1), ToolingUtils.java:222-282, CmdDiann.java:610-617).
 * Groups are getGroup2() (experiment_bioreplicate / exp_bioreplicate / experiment, InputLcmsFile.java:147-156)
 * in TreeMap order. GPF-DIA and DIA-Lib files are skipped, and so is every file that is not DIA or
 * DIA-Quant when there is DIA data (filterRun, ToolingUtils.java:284-292). Without any experiment or
 * bioreplicate, sample names come from the file names after their common prefix.
 *
 * @param entries Manifest rows from fragpipeManifestEntries
 * @return File content, header and one line per row, each ending in a newline
 */
def lfqExperimentAnnotation(entries) {
    def has_dia = entries.any { entry -> entry[3] in ['DIA', 'DIA-Quant'] }
    def skipped = { entry -> entry[3] in ['GPF-DIA', 'DIA-Lib'] || (has_dia && !(entry[3] in ['DIA', 'DIA-Quant'])) }
    def group2 = { entry ->
        if (entry[2] != null) {
            return entry[1].trim() ? "${entry[1]}_${entry[2]}".toString() : "exp_${entry[2]}".toString()
        }
        return entry[1]
    }
    def base_name = { path ->
        def parts = path.toString().tokenize('/')
        if (!parts) {
            error("lfqExperimentAnnotation: a manifest row has an empty LC-MS file path (row: ${entries.find { entry -> entry[0] == path }})")
        }
        def name = parts.last()
        name.lastIndexOf('.') >= 0 ? name.substring(0, name.lastIndexOf('.')) : name
    }
    def groups = entries.groupBy { entry -> group2.call(entry) }
    def lines = ['file\tsample\tsample_name\tcondition\treplicate']
    if (groups.size() == 1 && groups.keySet().first() == '') {
        def names = entries.collect { entry -> base_name.call(entry[0]) }.unique()
        // StringUtils.getCommonPrefix (commons-lang3)
        def prefix = names.inject(names[0]) { common, name ->
            def limit = Math.min(common.length(), name.length())
            def mismatch = (0..<limit).find { n -> common.charAt(n) != name.charAt(n) }
            common.substring(0, mismatch == null ? limit : mismatch)
        }
        entries.findAll { entry -> !skipped.call(entry) }.each { entry ->
            def sample_name = base_name.call(entry[0]).substring(prefix.length())
            def tt = sample_name.split('_')
            def condition = sample_name
            def replicate = 1
            if (tt.length == 3) {
                try {
                    replicate = Integer.parseInt(tt[2])
                    condition = "${tt[0]}_${tt[1]}"
                } catch (NumberFormatException _e) {
                    condition = "${tt[1]}_${tt[2]}"
                }
            } else if (tt.length == 2) {
                try {
                    replicate = Integer.parseInt(tt[1])
                    condition = tt[0]
                } catch (NumberFormatException _e) {
                    condition = tt[1]
                }
            }
            lines << [entry[0], entry[0], sample_name, condition, replicate].join('\t')
        }
    } else {
        groups.keySet().sort().each { group ->
            groups[group].findAll { entry -> !skipped.call(entry) }.each { entry ->
                def parts = entry[1].split('_')
                lines << [entry[0], entry[0], group, parts[0].trim(), entry[2] == null ? 1 : entry[2]].join('\t')
            }
        }
    }
    return lines.collect { line -> line + '\n' }.join('')
}

/**
 * LC-MS files SpecLibGen reads (CmdSpecLibGen.java:68-83, :172-185), in InputLcmsFile order: the group's
 * searched files, with MSFragger's <name>_uncalibrated.mzML for a .raw or timsTOF .d file (a DIA-Umpire /
 * diaTracer output takes the place of its source file). FragPipe always lists the uncalibrated mzML for
 * those, so its absence stops the run. Files outside .d/.mzML/.mzXML/.raw/.mgf stop the run as in FragPipe.
 *
 * @param searched List of [meta, spectra] searched by MSFragger (after DIA-Umpire / diaTracer)
 * @param uncalibrated Map of meta.id to MSFragger's _uncalibrated.mzML
 * @param files List of [meta, spectra] of the LC-MS table
 * @return Spectra files
 */
def speclibSpectraFiles(searched, uncalibrated, files) {
    def compatible = ['.d', '.mzml', '.mzxml', '.raw', '.mgf']
    def incompatible = files.findAll { item -> !compatible.any { ext -> item[1].name.toLowerCase().endsWith(ext) } }
    if (incompatible) {
        error("Spectral library generation with FragPipe-SpecLib is currently only compatible with ${compatible.join(', ')} input files. You can convert your data using msconvert program from ProteoWizard. Incompatible: ${incompatible.collect { item -> item[1].name }}")
    }
    def source = { id -> files.find { item -> item[0].id == id }?.getAt(1) }
    // FragPipe's path of a diaTracer output: the .d path with ".d" replaced by "_diatracer.mzML", next to the
    // .d (CmdDiaTracer.outputs, CmdDiaTracer.java:132-141), not the task dir that wrote it
    def diatracer_path = { f ->
        files.collect { item -> item[1].toAbsolutePath().normalize().toString() }
            .findAll { p -> p.toLowerCase().endsWith('.d') }
            .collect { p -> p.replace('.d', '_diatracer.mzML') }
            .find { p -> p.tokenize('/').last() == f.name }
    }
    // The group's LC-MS files in InputLcmsFile order (InputLcmsFile.compareTo, InputLcmsFile.java:116-119):
    // getGroup2(), data type, absolute normalised path as FragPipe holds it
    def sort_key = { item ->
        def experiment = fragpipeExperiment(item[0].experiment)
        def bioreplicate = item[0].bioreplicate
        def group2 = bioreplicate != null ? (experiment ? "${experiment}_${bioreplicate}" : "exp_${bioreplicate}") : experiment
        def table_path = source.call(item[0].id)
        def path = table_path ? table_path.toAbsolutePath().normalize().toString() : (diatracer_path.call(item[1]) ?: item[1].toAbsolutePath().normalize().toString())
        [group2.toString(), item[0].data_type ?: '', path]
    }
    return searched.sort(false) { a, b ->
        def ka = sort_key.call(a)
        def kb = sort_key.call(b)
        ka[0] <=> kb[0] ?: ka[1] <=> kb[1] ?: ka[2] <=> kb[2]
    }.collect { item ->
        def original = source.call(item[0].id) ?: item[1]
        def name = original.name.toLowerCase()
        if (!(name.endsWith('.raw') || name.endsWith('.d'))) {
            return item[1]
        }
        def spectra = uncalibrated[item[0].id]
        if (spectra == null) {
            error("SpecLibGen reads MSFragger's ${original.baseName}_uncalibrated.mzML for ${original.name} (CmdSpecLibGen.java:176-180), but MSFragger wrote none")
        }
        spectra
    }
}

/**
 * LC-MS files DIA-NN quantifies (CmdDiann.java:209-238): DIA and DIA-Quant files of the LC-MS table,
 * in table order. A DIA-Quant file keeps its path; without DIA-Umpire a DIA .raw file is replaced by
 * MSFragger's <name>_uncalibrated.mzML, whose absence stops the run.
 * Formats other than mzML/d/dia/raw stop the run (checkCompatibleFormats, CmdDiann.java:729-752).
 *
 * @param files List of [meta, spectra] of the samplesheet (meta.data_type set)
 * @param uncalibrated Map of meta.id to MSFragger's _uncalibrated.mzML
 * @param umpire_runs true when DIA-Umpire runs
 * @return DIA-NN input files ([] when there are none)
 */
def diannSpectraFiles(files, uncalibrated, umpire_runs) {
    def quantified = files.findAll { item -> item[0].data_type in ['DIA', 'DIA-Quant'] }
    def unsupported = quantified.collect { item -> item[1].name.toLowerCase().tokenize('.').last() }.unique()
        .findAll { ext -> !(ext in ['mzml', 'd', 'dia', 'raw']) }
    if (unsupported) {
        error("DIA-Quant can't work with '.${unsupported.join(', ')}' files. You can convert files using msconvert from ProteoWizard.")
    }
    return quantified.collect { item ->
        if (item[0].data_type == 'DIA-Quant' || umpire_runs || !item[1].name.toLowerCase().endsWith('.raw')) {
            return item[1]
        }
        def spectra = uncalibrated[item[0].id]
        if (spectra == null) {
            error("DIA-NN reads MSFragger's ${item[1].baseName}_uncalibrated.mzML for ${item[1].name} (CmdDiann.java:220), but MSFragger wrote none")
        }
        spectra
    }
}

/**
 * Which tables FPOP coADAPTr edits. FragPipe runs CmdFPOPcoadaptr after Report, IonQuant,
 * TMT-Integrator and DIA-NN (FragpipeRun.java:2487). With DIA-NN it edits DIA-NN's report.tsv
 * (FPOP_DiannReport_Editor); otherwise every psm.tsv of the work dir (FPOP_PSMwriter), which by
 * then carries the quantification columns of IonQuant (label-free), the TMT IonQuant isobaric
 * pass, or Philosopher report after FreeQuant/Labelquant; without quantification, the
 * Philosopher filter tables (CmdFPOPcoadaptr.java:429-433).
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @return 'diann', 'ionquant', 'isobaric', 'report' or 'filter'
 */
def fpopCoadaptrSource(json_text) {
    if (shouldRunTool(json_text, 'diann')) {
        return 'diann'
    }
    if (shouldRunTool(json_text, 'ionquant')) {
        return 'ionquant'
    }
    if (shouldRunTool(json_text, 'ionquant_isobaric')) {
        return 'isobaric'
    }
    if (shouldRunTool(json_text, 'freequant') || shouldRunTool(json_text, 'labelquant')) {
        return 'report'
    }
    return 'filter'
}

/**
 * Meta of a task whose tables FPOP coADAPTr replaces: meta.fpop_coadaptr names the source
 * (fpopCoadaptrSource), and that task's publishDir (conf/modules.config) then leaves the
 * unedited tables to FPOPCOADAPTR, which publishes the edited ones at the same paths.
 * FragPipe edits them in place, so only the edited version exists.
 *
 * @param meta Meta map of the task
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param sources The fpopCoadaptrSource values whose tables this task writes
 * @return meta, plus fpop_coadaptr: <source> when FPOP coADAPTr edits this task's tables
 */
def fpopCoadaptrMeta(meta, json_text, sources) {
    def source = fpopCoadaptrSource(json_text)
    return shouldRunTool(json_text, 'fpop_coadaptr') && source in sources ? meta + [fpop_coadaptr: source] : meta
}

/**
 * FPOP checks that need the LC-MS data types, raised before any task runs.
 * - FPOP TMT needs the second TMT-Integrator run, which FragPipe never runs on DIA data; FragPipe
 *   aborts with this message (FragpipeRun.java:2172-2193).
 * - Label-free FPOP reads IonQuant's combined_modified_peptide.tsv (CmdFpopQuant.java:62).
 *   FragPipe does not check this up front; its FPOP script then stops with "mod peptide file ...
 *   does not exist" (FragPipe_FPOP_Analysis.py:99-101). This pipeline check fails earlier.
 *
 * @param json_text JSON string from PARSEFRAGPIPEWORKFLOW tool_configs channel
 * @param data_types FragPipe data types of all LC-MS files (meta.data_type)
 * @return json_text unchanged
 */
def checkFpopTools(json_text, data_types) {
    def dia_types = ['DIA', 'DIA-Lib', 'DIA-Quant', 'GPF-DIA']
    if (shouldRunTool(json_text, 'tmtintegrator_unmod') && data_types.any { data_type -> data_type in dia_types }) {
        error("FPOP TMT Quant was requested on Downstream tab, but TMT-Integrator is not run. Please enable TMT-Integrator (FragPipe skips the second TMT-Integrator run for DIA data)")
    }
    if (shouldRunTool(json_text, 'fpop') && !shouldRunTool(json_text, 'tmtintegrator_unmod') && !shouldRunTool(json_text, 'ionquant')) {
        error("Label-free FPOP Quant reads IonQuant's combined_modified_peptide.tsv, but IonQuant does not run")
    }
    return json_text
}

/**
 * The work-dir fragpipe-files.fp-manifest FragPipe writes before a run (FragpipeRun.java:372 ->
 * TabWorkflow.manifestSave, TabWorkflow.java:1241-1250): one line per LC-MS file with the
 * InputLcmsFile experiment (sanitised, fragpipeExperiment), the bioreplicate and the data type,
 * written with the headless manifest writer (generateFragpipeManifest).
 * FragmentsExportMainClass (FragAnnotator-1.0.0.jar) reads it to find each group's psm.tsv:
 * <experiment>/ or <experiment>_<bioreplicate>/; when any experiment is blank (or SpecLibGen runs) it
 * reads only the work-dir root psm.tsv (ResultProcessor.goThroughFolder, inner_defined_empty_exp).
 * Spectra paths point to ../spectra/ beside the work dir. A .raw file keeps its name when MSFragger
 * wrote its <name>_uncalibrated.mzML / <name>_calibrated.mzML (the export reads those); a .raw
 * converted to mzML before the search is named as that <name>.mzML.
 *
 * @param files List of [meta, spectra] for every LC-MS file (meta.experiment, meta.bioreplicate, meta.data_type)
 * @param spectra_names File names of the staged spectra
 * @return Manifest text, one line per unique spectra file
 */
def matchedFragmentsManifest(files, spectra_names) {
    def entries = files.collect { meta, spectra ->
        def name = spectra.name
        def msfragger_mzml = ["${spectra.baseName}_uncalibrated.mzML", "${spectra.baseName}_calibrated.mzML"].collect { n -> n.toString() }
        if (name.toLowerCase().endsWith('.raw') && !msfragger_mzml.any { n -> n in spectra_names }) {
            name = "${spectra.baseName}.mzML"
        }
        ["../spectra/${name}".toString(), fragpipeExperiment(meta.experiment), meta.bioreplicate, meta.data_type]
    }.unique().sort { entry -> entry[0] }
    return generateFragpipeManifest(entries)
}

/**
 * MBG --nopasef value: TabWorkflow.getInputDataType() != ImMsTimsTof (CmdMBGMatch.java:196-198).
 * Any .d file means timsTOF; otherwise any .raw file means regular MS; otherwise the workflow's
 * IM-MS radio button decides (TabWorkflow.java:521-538).
 *
 * @param file_names Names of every LC-MS file in the run, as the FragPipe file table lists them
 * @param im_ms_radio workflow.input.data-type.im-ms
 * @return "true" or "false"
 */
def mbgNoPasef(file_names, im_ms_radio) {
    def names = file_names.collect { name -> name.toString().toLowerCase() }
    def timstof = names.any { name -> name.endsWith('.d') } ? true
        : names.any { name -> name.endsWith('.raw') } ? false
        : im_ms_radio
    return timstof ? 'false' : 'true'
}

/**
 * Spectra file O-Pair's filelist names for one LC-MS file (CmdOPair.java:107-114): the
 * MSFragger <base>_calibrated.mzML when calibrated mzML is written, else <base>_uncalibrated.mzML
 * for .raw/.d, else the input file itself.
 *
 * @param spectra_kind tool_configs opair.spectra ('calibrated' or 'uncalibrated_or_input')
 * @param original_name Name of the LC-MS file in the FragPipe file table
 * @return 'calibrated', 'uncalibrated' or 'input'
 */
def opairSpectraSource(spectra_kind, original_name) {
    if (spectra_kind == 'calibrated') {
        return 'calibrated'
    }
    def name = original_name.toString().toLowerCase()
    return (name.endsWith('.raw') || name.endsWith('.d')) ? 'uncalibrated' : 'input'
}

/**
 * FragPipe 24.0's templates PARSEFRAGPIPEWORKFLOW stages: fragger_closed.params (the MSFragger tab's
 * closed-search defaults) and diaumpire_se.params (the DIA-Umpire template), kept in the module's assets/.
 *
 * @param project_dir Nextflow projectDir path
 * @return [fragger_closed.params, diaumpire_se.params] files; a missing file stops the run
 */
def parseFragpipeWorkflowTemplates(project_dir) {
    return ['fragger_closed.params', 'diaumpire_se.params'].collect { name ->
        file("${project_dir}/modules/local/parsefragpipeworkflow/assets/${name}", checkIfExists: true)
    }
}

/**
 * The qpx export (.h5mu / .h5ad / QPX Parquet) a FragPipe run gets, checked before any task runs.
 *
 * - 'tmt': TMT-Integrator runs; its abundance_*.tsv reports go to TMTANNDATA.
 * - 'lfq': label-free IonQuant runs (ionquant.run-ionquant=true, no TMT-Integrator); its
 *   combined_protein.tsv / combined_ion.tsv go to QPXFRAGPIPE, which needs the MaxLFQ intensities.
 *   IonQuant writes them only with at least 3 experiment+bioreplicate combinations ("Add MaxLFQ",
 *   QuantPanelLabelfree.java:316), i.e. FragPipe groups, so fewer stop the run here instead of after the search.
 * - '': no export (--skip_qpx_export, the TMT label check, or a workflow without IonQuant or TMT-Integrator).
 *
 * @param skip           params.skip_qpx_export
 * @param tmt_mode       params.tmt_mode ('labelcheck' has no export)
 * @param workflow_props FragPipe workflow properties (fragpipeWorkflowProperties)
 * @param spectra_files  List of [meta, spectra] of every LC-MS file (meta.group = FragPipe group)
 * @return 'lfq', 'tmt' or ''
 */
def qpxFragpipeExport(skip, tmt_mode, workflow_props, spectra_files) {
    if (skip || tmt_mode == 'labelcheck') {
        return ''
    }
    if (workflow_props.getProperty('tmtintegrator.run-tmtintegrator') == 'true') {
        return 'tmt'
    }
    if (workflow_props.getProperty('ionquant.run-ionquant') != 'true') {
        return ''
    }
    def min_groups = 3
    def groups = spectra_files.collect { item -> item[0].group }.unique().sort()
    if (groups.size() < min_groups) {
        error("qpx export: the samplesheet has ${groups.size()} experiment+bioreplicate combination(s) ${groups}, fewer than the ${min_groups} IonQuant needs to write the MaxLFQ intensities QPXFRAGPIPE converts. Add experiments, or rerun with --skip_qpx_export true.")
    }
    return 'lfq'
}

/**
 * TMT-Integrator's log2transformed setting from the tmt-integrator-conf.yml text TMT-Integrator reads
 * ("  log2transformed: <value>" under "tmtintegrator:", TmtiConfig.java:332-338): the tool configs'
 * tmtintegrator args in pipeline mode, FragPipe's <workdir>/tmt-integrator-conf.yml in headless mode
 * (CmdTmtIntegrator.java:160).
 * Returns null when the text has no log2transformed line; the caller raises the error (an error raised in
 * an included function called from an operator closure reaches the console as InvocationTargetException).
 *
 * @param config_text tmt-integrator-conf.yml text
 * @return 'true', 'false', or null when the setting is absent
 */
def tmtLog2Transformed(config_text) {
    def line = config_text.toString().readLines().find { l -> l ==~ /\s*log2transformed:.*/ }
    return line == null ? null : (line.split(':', 2)[1].trim().toLowerCase() == 'true' ? 'true' : 'false')
}

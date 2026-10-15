//
// Subworkflow with functionality specific to the nf-core/msproteomics pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/msproteomics ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/msproteomics/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --mode <diann|fragpipe> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    emit:
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def tools = params.mode == 'diann'
        ? ["DIA-NN (Demichev et al. 2020),", "quantms-utils (Dai et al. 2024),", "MSstats (Kohler et al. 2023),", "iq (Pham et al. 2020),", "pmultiqc,"]
        : ["FragPipe (Yu et al. 2023): MSFragger (Kong et al. 2017), MSBooster (Yang et al. 2023), Percolator (Kall et al. 2007), Philosopher (da Veiga Leprevost et al. 2020), IonQuant (Yu et al. 2021) and TMT-Integrator (Djomehri et al. 2020),"]
    return (["Tools used in the workflow included:"] + tools + ["MultiQC (Ewels et al. 2016)", "."]).join(' ').trim()
}

def toolBibliographyText() {
    def references = params.mode == 'diann'
        ? [
            "<li>Demichev, V., Messner, C. B., Vernardis, S. I., Lilley, K. S., & Ralser, M. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nature Methods, 17(1), 41-44. doi: 10.1038/s41592-019-0638-x</li>",
            "<li>Kohler, D., Staniak, M., Tsai, T. H., Huang, T., Shulman, N., Bernhardt, O. M., MacLean, B. X., Nesvizhskii, A. I., Reiter, L., Sabido, E., Choi, M., & Vitek, O. (2023). MSstats Version 4.0: Statistical Analyses of Quantitative Mass Spectrometry-Based Proteomic Experiments with Chromatography-Based Quantification at Scale. Journal of Proteome Research, 22(5), 1466-1482. doi: 10.1021/acs.jproteome.2c00834</li>",
            "<li>Pham, T. V., Henneman, A. A., & Jimenez, C. R. (2020). iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics. Bioinformatics, 36(8), 2611-2613. doi: 10.1093/bioinformatics/btz961</li>"
        ]
        : [
            "<li>Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., & Nesvizhskii, A. I. (2017). MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics. Nature Methods, 14(5), 513-520. doi: 10.1038/nmeth.4256</li>",
            "<li>da Veiga Leprevost, F., Haynes, S. E., Avtonomov, D. M., Chang, H. Y., Shanmugam, A. K., Mellacheruvu, D., Kong, A. T., & Nesvizhskii, A. I. (2020). Philosopher: a versatile toolkit for shotgun proteomics data analysis. Nature Methods, 17(9), 869-870. doi: 10.1038/s41592-020-0912-y</li>",
            "<li>Yu, F., Haynes, S. E., & Nesvizhskii, A. I. (2021). IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Molecular & Cellular Proteomics, 20, 100077. doi: 10.1016/j.mcpro.2021.100077</li>"
        ]
    references << "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
    return references.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

/**
 * Convert an organism name to standard binomial nomenclature.
 *
 * @param organism Raw organism name (e.g., 'human', 'Homo_sapiens', 'hs')
 * @return Standard organism name (e.g., 'Homo sapiens'); unknown names are returned unchanged
 */
def convertOrganismToStandardName(organism) {
    def org = organism.toLowerCase().strip().replaceAll("\\s+", "_")
    if (org == "homo_sapiens" || org == "human" || org == "hs") {
        return "Homo sapiens"
    } else if (org == "mus_musculus" || org == "mouse" || org == "mm") {
        return "Mus musculus"
    } else if (org == "saccharomyces_cerevisiae" || org == "yeast" || org == "sc") {
        return "Saccharomyces cerevisiae"
    }
    return organism
}

/**
 * The iProteome entry selected with --proteome (params.iproteomes, conf/iproteomes.config).
 *
 * @param pipeline_params Pipeline params map
 * @return The entry map, or null when --proteome is not set
 */
def resolveProteome(pipeline_params) {
    def key = pipeline_params.proteome
    if (!key) {
        return null
    }
    if (pipeline_params.iproteomes_ignore) {
        error("--proteome '${key}' needs the iProteomes config: remove --iproteomes_ignore")
    }
    if (!pipeline_params.iproteomes_base) {
        error("--proteome '${key}' needs --iproteomes_base: the root directory of the iProteomes store")
    }
    def entry = pipeline_params.iproteomes ? pipeline_params.iproteomes[key] : null
    if (!entry) {
        error("--proteome '${key}' is not an iProteome; available: ${pipeline_params.iproteomes?.keySet()?.sort()?.join(', ')}")
    }
    if (!entry.fasta) {
        error("iProteome '${key}' declares no fasta")
    }
    // Entry paths are relative to --iproteomes_base, joined here so that -c files can set the base
    def base = pipeline_params.iproteomes_base.toString().replaceAll('/+$', '')
    return entry + [key: key] + ['fasta', 'fragpipe_fasta', 'diann_speclib']
        .findAll { attribute -> entry[attribute]?.toString()?.trim() }
        .collectEntries { attribute -> [(attribute): "${base}/${entry[attribute].toString().trim()}".toString()] }
}

/**
 * Hosted index (iProteome attribute 'diann_speclib' or 'fragpipe_fasta') of an entry.
 * An entry that declares an index whose path is empty, missing or unreadable is an error;
 * an entry that does not declare it has none to use.
 *
 * @param entry iProteome entry map, or null
 * @param attribute Index attribute name
 * @param key iProteome key (for messages)
 * @return Index path, or null when the entry does not declare the index
 */
def proteomeIndex(entry, attribute, key) {
    if (!entry || !entry.containsKey(attribute)) {
        return null
    }
    def path = entry[attribute]?.toString()?.trim()
    if (!path) {
        error("iProteome '${key}' declares '${attribute}' but its path is empty")
    }
    def index = file(path)
    if (!index.exists() || !index.canRead()) {
        error("iProteome '${key}' declares '${attribute}' but ${path} is missing or unreadable")
    }
    return path
}

/**
 * The workflow's database.db-path against the FASTA the run searches (--database or the iProteome).
 *
 * FragPipe searches the workflow's own FASTA: the Database tab field "db-path" (TabDatabase.java:142) is
 * the path every tool gets (fastaPath = checkFasta(tabRun, configDb), configureTaskGraph(..., fastaPath, ...),
 * FragpipeRun.java:327, :344). The pipeline takes --database instead: the workflow's db-path is a path on the
 * machine that saved it (not readable from Tower/Batch), and renamed copies of one FASTA are normal.
 *   - db-path = FragPipe's databases folder, as the finalizer writes it into the second-pass workflow
 *     (FragpipeRun.java:671-672, FragpipeLocations.getDirDatabases() = <install>/databases): FragPipe stops,
 *     a folder is no regular file (TabDatabase.java:274-276 posts an empty NoteConfigDatabase, then
 *     checkFasta logs "FASTA file path is empty or the file is corrupted", FragpipeRun.java:1046-1049).
 *     The run stops the same way.
 *   - another file name: a warning naming both files; --database is searched.
 *   - blank db-path: --database is searched.
 *
 * @param workflow_props Properties of the .workflow file
 * @param databases      FASTA paths or URLs the run searches (--database / iProteome, hosted FragPipe FASTA)
 */
def checkWorkflowDatabase(workflow_props, databases) {
    def db_path = (workflow_props.getProperty('database.db-path') ?: '').trim()
    if (!db_path) {
        return
    }
    // the folder's last path segment is "databases", whatever the install (FragPipe-24.0\databases on Windows, /opt/fragpipe/databases)
    if (db_path.replace('\\', '/').replaceAll(/\/+$/, '').tokenize('/').last().equalsIgnoreCase('databases')) {
        error("FASTA file path is empty or the file is corrupted (database.db-path=${db_path} is FragPipe's databases folder, as in a second-pass workflow; set it to the FASTA file). FragpipeRun.java:1046-1049, TabDatabase.java:274-276")
    }
    def name = db_path.tokenize('/\\').last()
    def names = databases.collect { db -> db.toString().tokenize('/').last() }.unique()
    if (!(name in names)) {
        log.warn("The workflow's database.db-path (${db_path}) names another file than the FASTA the run searches (${names.join(', ')}); --database replaces the workflow's db-path")
    }
}

/**
 * The protein FASTA searched. Priority: --database > the iProteome fasta (--proteome) > error.
 *
 * @param entry iProteome entry map, or null
 * @param pipeline_params Pipeline params map
 * @return FASTA path or URL
 */
def resolveDatabase(entry, pipeline_params) {
    def database = pipeline_params.database ?: entry?.fasta
    if (!database) {
        error("Needs a protein FASTA: set --database or --proteome (iProteomes, --iproteomes_base)")
    }
    return database
}

/**
 * Organism for the SDRF bookkeeping: --organism > the iProteome organism > '' (empty: the SDRF
 * writer guesses it from the searched FASTA's OS= headers with FragPipe's rule, the same
 * implementation as FragPipe's sdrf.tsv, bin/fragpipe_sdrf.py guess_organism).
 *
 * @param entry iProteome entry map, or null
 * @param pipeline_params Pipeline params map
 * @return Standard organism name, or '' to guess from the database
 */
def resolveOrganism(entry, pipeline_params) {
    def organism = pipeline_params.organism ?: entry?.organism
    return organism ? convertOrganismToStandardName(organism) : ''
}

/**
 * Normalised DIA-NN library setting: numbers compare as numbers (18 == 18.0), comma-separated
 * modification lists as sets (order-insensitive), anything else as trimmed text.
 *
 * @param value Setting value (null = unset)
 * @return Normalised text; '' when unset
 */
def normaliseLibrarySetting(value) {
    def text = value == null ? '' : value.toString().trim()
    if (text ==~ /-?\d+(\.\d+)?/) {
        return new BigDecimal(text).stripTrailingZeros().toPlainString()
    }
    return text.contains(',') ? text.split(',').collect { item -> item.trim() }.sort().join(',') : text
}

/**
 * DIA-NN library settings of a run (params names) that decide the content of an in-silico library.
 *
 * @param pipeline_params Pipeline params map
 * @return Normalised settings map
 */
def diannLibrarySettings(pipeline_params) {
    return [
        'max_mods', 'allowed_missed_cleavages',
        'min_peptide_length', 'max_peptide_length', 'min_precursor_charge', 'max_precursor_charge',
        'diann_min_pr_mz', 'diann_max_pr_mz', 'diann_min_fr_mz', 'diann_max_fr_mz',
        'diann_library_mass_acc', 'diann_library_ms1_acc', 'additional_insilico_args'
    ].collectEntries { name -> [(name): normaliseLibrarySetting(pipeline_params[name])] } +
        [
            enzyme: normaliseLibrarySetting(diannEnzymeCut(pipeline_params.enzyme ?: 'Trypsin/P')),
            fixed_mods: diannModFlags(pipeline_params.fixed_mods),
            variable_mods: diannModFlags(pipeline_params.variable_mods),
            met_excision: pipeline_params.met_excision ? 'true' : 'false'
        ]
}

/**
 * DIA-NN --cut rule of an --enzyme value. Port of quantms-utils 0.0.23 (the version of
 * modules/nf-core/quantmsutils/dianncfg/environment.yml) quantmsutils/diann/dianncfg.py,
 * _ENZYME_SPECIFICITY (lines 126-133) and enzyme_cut (line 136).
 *
 * @param enzyme --enzyme value
 * @return DIA-NN --cut rule
 */
def diannEnzymeCut(enzyme) {
    def specificity = [
        'Trypsin': 'K*,R*,!*P', 'Trypsin/P': 'K*,R*', 'Arg-C': 'R*,!*P',
        'Asp-N': '*B,*D', 'Chymotrypsin': 'F*,W*,Y*,L*,!*P', 'Lys-C': 'K*,!*P'
    ]
    if (!specificity.containsKey(enzyme)) {
        error("enzyme '${enzyme}' has no DIA-NN cut rule in quantms-utils dianncfg (known: ${specificity.keySet()})")
    }
    return specificity[enzyme]
}

/**
 * Monoisotopic delta mass of the Unimod modifications dianncfg can look up by name:
 * <delta mono_mass> of the first <mod title=...> in sdrf_pipelines 0.0.32 (the version of
 * modules/nf-core/quantmsutils/dianncfg/environment.yml) sdrf_pipelines/openms/unimod.xml,
 * read by UnimodDatabase (openms/unimod.py:56-110) as dianncfg.py:60-62 uses it.
 */
def diannUnimodMasses() {
    return [
        'Acetyl': '42.010565', 'Carbamidomethyl': '57.021464', 'Carbamyl': '43.005814',
        'Deamidated': '0.984016', 'Phospho': '79.966331', 'Methyl': '14.01565',
        'Oxidation': '15.994915', 'Dimethyl': '28.0313', 'Trimethyl': '42.04695'
    ]
}

/**
 * Normalised "<mass>@<site>" set of a --fixed_mods / --variable_mods value, converted as
 * quantms-utils 0.0.23 dianncfg.py convert_mod (lines 51-123) builds the DIA-NN modification
 * flags 'name,mass,site': name = text before the first space (line 61/97), site = text in the
 * first parentheses of the rest (line 71/107), 'Protein N-term' -> '*n', 'N-term' -> 'n' (lines 72-75).
 *
 * @param mods Comma-separated 'Name (Site)' list, or null
 * @return Sorted ';'-joined '<mass>@<site>' set; '' when none
 */
def diannModFlags(mods) {
    def text = mods == null ? '' : mods.toString()
    if (text == '') {
        return ''
    }
    return text.split(',').collect { mod ->
        def name = mod.split(' ')[0]
        def mass = diannUnimodMasses()[name]
        if (!mass) {
            error("modification '${mod}': '${name}' has no recorded Unimod mass (known: ${diannUnimodMasses().keySet()})")
        }
        def site = (mod.split(' ', 2)[1] =~ /\((.*?)\)/)[0][1]
        site = site == 'Protein N-term' ? '*n' : (site == 'N-term' ? 'n' : site)
        "${normaliseLibrarySetting(mass)}@${site}".toString()
    }.sort().join(';')
}

/**
 * Settings an iProteome's DIA-NN library was built with, read from the diann_library_settings
 * object of its meta.json (<iproteomes_base>/<key>/meta.json), expressed as the DIA-NN flags (--cut, mass@site) the run's params convert to.
 *
 * @param entry iProteome entry map from resolveProteome
 * @param pipeline_params Pipeline params map
 * @return Normalised settings map, same keys as diannLibrarySettings
 */
def readDiannLibraryRecord(entry, pipeline_params) {
    def meta_file = file("${pipeline_params.iproteomes_base}/${entry.key}/meta.json")
    if (!meta_file.exists()) {
        error("iProteome '${entry.key}' declares 'diann_speclib' but ${meta_file} is missing: the settings the library was built with are unknown")
    }
    def meta = new groovy.json.JsonSlurper().parseText(meta_file.text)
    def lib = meta.diann_library_settings
    if (!lib) {
        error("${meta_file} has no diann_library_settings object: the settings the library was built with are unknown")
    }
    def modFlags = { mods ->
        mods.collect { mod -> "${normaliseLibrarySetting(mod.mass)}@${mod.residues}".toString() }.sort().join(';')
    }
    return [
        enzyme: lib.cut,
        fixed_mods: modFlags.call(lib.fixed_mods),
        variable_mods: modFlags.call(lib.var_mods),
        max_mods: lib.max_var_mods,
        allowed_missed_cleavages: lib.missed_cleavages,
        min_peptide_length: lib.min_pep_len, max_peptide_length: lib.max_pep_len,
        min_precursor_charge: lib.min_pr_charge, max_precursor_charge: lib.max_pr_charge,
        diann_min_pr_mz: lib.min_pr_mz, diann_max_pr_mz: lib.max_pr_mz,
        diann_min_fr_mz: lib.min_fr_mz, diann_max_fr_mz: lib.max_fr_mz,
        diann_library_mass_acc: lib.mass_acc, diann_library_ms1_acc: lib.mass_acc_ms1,
        met_excision: lib.met_excision ? 'true' : 'false',
        additional_insilico_args: ''
    ].collectEntries { name, value -> [(name): normaliseLibrarySetting(value)] }
}

/**
 * The DIA-NN spectral library. Priority: --diann_speclib > the iProteome library (--proteome,
 * no --database: the library belongs to the iProteome FASTA) > null (built in the pipeline).
 * The hosted library is used only when the run's digest/modification settings equal the
 * settings recorded in the iProteome's meta.json.
 *
 * @param entry iProteome entry map, or null
 * @param pipeline_params Pipeline params map
 * @return Library path, or null when the pipeline builds it
 */
def resolveDiannSpeclib(entry, pipeline_params) {
    if (pipeline_params.diann_speclib) {
        return pipeline_params.diann_speclib
    }
    if (pipeline_params.database) {
        return null
    }
    def library = proteomeIndex(entry, 'diann_speclib', pipeline_params.proteome)
    if (library) {
        def recorded = readDiannLibraryRecord(entry, pipeline_params)
        def run = diannLibrarySettings(pipeline_params)
        def differing = recorded.keySet().findAll { name -> recorded[name] != run[name] }
        if (differing) {
            error("iProteome '${pipeline_params.proteome}' DIA-NN library was built with different settings than this run: " +
                differing.collect { name -> "${name} library='${recorded[name]}' run='${run[name]}'" }.join('; ') +
                ". Set the run to the library's settings, or pass --diann_speclib / --database to use another library or have the pipeline build it.")
        }
    }
    return library
}

/**
 * The protein FASTA FragPipe searches (database.db-path). Priority: --database > the iProteome
 * fragpipe_fasta (decoys and contaminants added; not with --add_philosopher_decoy /
 * --add_philosopher_contaminants, which update the plain FASTA) > the iProteome fasta > error.
 * FragPipe headless searches database.db-path as given: a FASTA without decoys stops in
 * FragPipe's decoy check (FragpipeRun.checkDbConfig).
 *
 * @param entry iProteome entry map, or null
 * @param pipeline_params Pipeline params map
 * @return FASTA path or URL
 */
def resolveFragPipeDatabase(entry, pipeline_params) {
    def update = pipeline_params.add_philosopher_decoy || pipeline_params.add_philosopher_contaminants
    if (!pipeline_params.database && !update) {
        def hosted = proteomeIndex(entry, 'fragpipe_fasta', pipeline_params.proteome)
        if (hosted) {
            return hosted
        }
    }
    return resolveDatabase(entry, pipeline_params)
}

/**
 * FragPipe database channels. The database is searched as given (FragpipeRun.checkDbConfig only
 * checks its decoys) unless --add_philosopher_decoy and/or --add_philosopher_contaminants ask for
 * FragPipe's database update (CmdDatabaseUpdate) first.
 *
 * @param entry iProteome entry map, or null
 * @param pipeline_params Pipeline params map
 * @return Map with 'ch_fasta' (FASTA to update; empty unless a database update is asked for) and
 *         'ch_prebuilt_db' (FASTA searched as given; empty when ch_fasta is used)
 */
def resolveFragPipeDatabaseChannels(entry, pipeline_params) {
    def ch_database = channel
        // glob: false keeps the query string of a UniProt REST URL (nf-core/msproteomics e5ce83c)
        .fromPath(resolveFragPipeDatabase(entry, pipeline_params), glob: false, checkIfExists: true)
        .map { db -> [[id: 'database'], db] }
    return (pipeline_params.add_philosopher_decoy || pipeline_params.add_philosopher_contaminants)
        ? [ch_fasta: ch_database, ch_prebuilt_db: channel.empty()]
        : [ch_fasta: channel.empty(), ch_prebuilt_db: ch_database]
}

/**
 * Generate the SDRF bookkeeping JSON from pipeline params.
 *
 * @param pipeline_params Pipeline params map
 * @param organism Standard organism name (resolveOrganism)
 * @return JSON string with SDRF metadata
 */
def generateSdrfBookkeeping(pipeline_params, organism) {
    def params_map = [
        organism:                  organism,
        enzyme:                    pipeline_params.enzyme ?: 'Trypsin/P',
        fixed_mods:                pipeline_params.fixed_mods ?: '',
        variable_mods:             pipeline_params.variable_mods ?: '',
        instrument:                '',
        acquisition_method:        pipeline_params.mode == 'diann' ? 'DIA' : 'DDA',
        precursor_mass_tolerance:  "${pipeline_params.precursor_mass_tolerance} ${pipeline_params.precursor_mass_tolerance_unit}",
        fragment_mass_tolerance:   "${pipeline_params.fragment_mass_tolerance} ${pipeline_params.fragment_mass_tolerance_unit}",
        dissociation_method:       pipeline_params.dissociation_method ?: 'HCD'
    ]
    return groovy.json.JsonOutput.toJson(params_map)
}

/**
 * Get the bundled TMT label check workflow file for a TMT type.
 *
 * @param tmt_type TMT type string (e.g., 'TMT18', 'TMT10')
 * @param project_dir Nextflow projectDir path
 * @return Path string of assets/TMT-labelcheck-<mass>.workflow
 */
def getTmtLabelcheckWorkflow(tmt_type, project_dir) {
    if (!tmt_type) {
        error("--tmt_type is required for TMT workflows (e.g., --tmt_type TMT18)")
    }
    def mass_map = [
        'TMT6': 229, 'TMT10': 229, 'TMT11': 229,
        'TMT16': 304, 'TMT18': 304, 'TMT35': 304, 'TMTPRO': 304
    ]
    def mass = mass_map[tmt_type.toUpperCase()]
    if (!mass) {
        error("Unknown tmt_type '${tmt_type}'. Valid values: ${mass_map.keySet().join(', ')}")
    }
    return "${project_dir}/assets/TMT-labelcheck-${mass}.workflow"
}

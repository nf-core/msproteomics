/*
 * WorkflowUtils: Pure Groovy utility functions for the msproteomics pipeline.
 *
 * These functions handle input parsing, organism name resolution, database
 * resolution, workflow validation, and FragPipe manifest generation.
 * They are auto-loaded from lib/ and can be called as
 * WorkflowUtils.methodName() from any Nextflow script.
 */
class WorkflowUtils {

    /**
     * Get a single unique value from a specified column in a delimited file.
     *
     * @param filePath Path to the input file
     * @param separator Field separator (e.g., '\t', ',')
     * @param columnName Name of the column to extract
     * @param defaultValue Optional default if column is empty (default: throws error)
     * @return Single unique value from the specified column
     */
    static String getUniqueColumnValue(filePath, String separator, String columnName, String defaultValue = '__NO_DEFAULT__') {
        def lines = new File(filePath.toString()).readLines()
        def header = lines[0].split(separator).collect { it.trim() }
        def col_idx = header.findIndexOf { it.trim().toLowerCase() == columnName.trim().toLowerCase() }

        if (col_idx == -1) {
            throw new RuntimeException("No column named '${columnName}' found in the input file")
        }

        def unique_values = lines.tail().collect { it.split(separator)[col_idx]?.trim() ?: '' }.unique()
        def non_empty_values = unique_values.findAll { it != '' }
        if (non_empty_values.size() > 1) {
            throw new RuntimeException("Multiple values found in column '${columnName}': ${non_empty_values}")
        }
        def ov = non_empty_values.size() > 0 ? non_empty_values[0] : ''
        if (ov == '') {
            if (defaultValue != '__NO_DEFAULT__') {
                return defaultValue
            } else {
                throw new RuntimeException("No unique value found in column '${columnName}'")
            }
        }
        return ov
    }

    /**
     * Convert organism name to standard binomial nomenclature.
     *
     * @param organism Raw organism name (e.g., 'human', 'Homo_sapiens', 'hs')
     * @return Standard organism name (e.g., 'Homo sapiens')
     */
    static String convertOrganismToStandardName(String organism) {
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
     * Resolve FragPipe database path and type.
     * Priority: philosopher_database > database > databases[organism].philosopher_database > databases[organism].database > error
     *
     * @param organism Organism name (e.g., 'Homo sapiens')
     * @param params Pipeline params object
     * @return Map with 'path' (String) and 'is_philosopher' (boolean) keys
     */
    static Map resolveFragPipeDatabase(String organism, params) {
        if (params.philosopher_database) {
            return [path: params.philosopher_database, is_philosopher: true]
        }
        if (params.database) {
            return [path: params.database, is_philosopher: false]
        }
        if (organism && params.databases && params.databases[organism]?.philosopher_database) {
            return [path: params.databases[organism].philosopher_database, is_philosopher: true]
        }
        if (organism && params.databases && params.databases[organism]?.database) {
            return [path: params.databases[organism].database, is_philosopher: false]
        }
        throw new RuntimeException(
            "No database provided. Use --philosopher_database, --database, " +
            "or ensure organism '${organism}' has a configured database."
        )
    }

    /**
     * Resolve FragPipe database and create appropriate Nextflow channels.
     *
     * @param organism Organism name
     * @param params Pipeline params object
     * @return Map with 'ch_fasta', 'ch_prebuilt_db' channels and 'resolved_db' info
     */
    static Map resolveFragPipeDatabaseChannels(String organism, params) {
        def resolved_db = resolveFragPipeDatabase(organism, params)
        def ch_fasta
        def ch_prebuilt_db

        if (resolved_db.is_philosopher) {
            ch_fasta = nextflow.Channel.empty()
            ch_prebuilt_db = nextflow.Channel
                .fromPath(resolved_db.path, checkIfExists: true)
                .map { db -> [[id: "${organism}_philosopher_database"], db] }
        } else {
            ch_fasta = nextflow.Channel
                .fromPath(resolved_db.path, checkIfExists: true)
                .map { db -> [[id: 'database'], db] }
            ch_prebuilt_db = nextflow.Channel.empty()
        }

        return [ch_fasta: ch_fasta, ch_prebuilt_db: ch_prebuilt_db, resolved_db: resolved_db]
    }

    /**
     * Parse SDRF channel to get sample input files with metadata.
     * Returns a channel of [meta, file_path] tuples.
     *
     * @param ch_sdrf Channel containing path to SDRF TSV file
     * @param sampleDir Optional sample directory override
     * @return Channel of [meta, file_path] tuples
     */
    static parseSdrfInputsChannel(ch_sdrf, sampleDir) {
        return ch_sdrf
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def file_name = row['comment[data file]']
                def file_uri = row['comment[file uri]'] ?: file_name
                def raw_source_name = row['source name']
                def sample_name = row['comment[sample name]']
                def source_name = (raw_source_name && raw_source_name != 'not available') ?
                    raw_source_name :
                    (sample_name ?: file_name.replaceAll(/\.(raw|mzML|mzml|RAW|d)$/, ''))

                def file_path = sampleDir ?
                    nextflow.Nextflow.file("${sampleDir}/${file_name}") :
                    nextflow.Nextflow.file(file_uri)

                def fraction_id = row['comment[fraction identifier]']
                def is_fractionated = fraction_id && fraction_id != 'not available'
                def meta = [
                    id:    is_fractionated ? "${source_name}_f${fraction_id}" : source_name,
                    group: source_name
                ]
                [ meta, file_path ]
            }
            .filter { meta, file_path ->
                if (!file_path.exists()) {
                    log.warn "File not found: ${file_path}, skipping..."
                    return false
                }
                return true
            }
    }

    /**
     * Generate SDRF bookkeeping JSON from pipeline params.
     *
     * @param params Pipeline params object
     * @return JSON string with SDRF metadata
     */
    static String generateSdrfBookkeeping(params) {
        def params_map = [
            organism:                  params.organism ?: 'Homo sapiens',
            enzyme:                    params.enzyme ?: 'Trypsin/P',
            fixed_mods:                params.fixed_mods ?: '',
            variable_mods:             params.variable_mods ?: '',
            instrument:                '',
            acquisition_method:        params.mode == 'diann' ? 'DIA' : 'DDA',
            precursor_mass_tolerance:  "${params.precursor_mass_tolerance ?: 20} ppm",
            fragment_mass_tolerance:   "${params.fragment_mass_tolerance ?: 20} ppm",
            dissociation_method:       params.dissociation_method ?: 'HCD'
        ]
        return groovy.json.JsonOutput.toJson(params_map)
    }

    /**
     * Get the TMT label check workflow file path based on TMT type.
     *
     * @param tmt_type TMT type string (e.g., 'TMT18', 'TMT10')
     * @param projectDir Nextflow projectDir path
     * @return Path to the TMT labelcheck workflow file
     */
    static String getTmtLabelcheckWorkflow(String tmt_type, projectDir) {
        if (!tmt_type) {
            throw new RuntimeException("ERROR: --tmt_type is required for TMT workflows (e.g., --tmt_type TMT18)")
        }
        def mass_map = [
            'TMT6': 229, 'TMT10': 229, 'TMT11': 229,
            'TMT16': 304, 'TMT18': 304, 'TMTPRO': 304
        ]
        def mass = mass_map[tmt_type.toUpperCase()]
        if (!mass) {
            throw new RuntimeException("ERROR: Unknown tmt_type '${tmt_type}'. Valid values: ${mass_map.keySet().join(', ')}")
        }
        return "${projectDir}/assets/TMT-labelcheck-${mass}.workflow"
    }

    /**
     * Parse TMT settings and data type from a FragPipe .workflow file.
     *
     * @param workflowFilePath Path to the FragPipe .workflow file
     * @return Map with 'is_tmt', 'plex_type', and 'data_type' keys
     */
    static Map parseWorkflowFileTmtSettings(workflowFilePath) {
        def text = new File(workflowFilePath.toString()).text
        def is_tmt = text.contains('tmtintegrator.run-tmtintegrator=true')
        def match = (text =~ /tmtintegrator\.channel_num=(.+)/)
        def plex_type = match ? match[0][1].trim() : null
        def data_type = text.contains('workflow.input.data-type.im-ms=true') ? 'DDA+' : 'DDA'
        return [is_tmt: is_tmt, plex_type: plex_type, data_type: data_type]
    }

    /**
     * Generate a FragPipe manifest from a list of [group, filename] entries.
     *
     * @param entries List of [group, filename] pairs
     * @param dataType Data type string (e.g., 'DDA', 'DDA+')
     * @return Tab-separated manifest string
     */
    static String generateFragpipeManifest(entries, String dataType) {
        def group_counts = entries.countBy { it[0] }
        return entries.collect { entry ->
            def group = entry[0]
            def fname = entry[1]
            def biorep = group_counts[group] > 1 ? '' : '1'
            "${fname}\t${group}\t${biorep}\t${dataType}"
        }.join('\n')
    }
}

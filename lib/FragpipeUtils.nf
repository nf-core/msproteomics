/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FragPipe Utility Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared helper functions for FragPipe Nextflow subworkflows.
    Import via: include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

    These functions accept a JSON text string (from tool_configs value channel)
    and parse it on demand. JSON parsing of the ~2-5KB config is microseconds,
    so the overhead of repeated parsing is negligible.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.json.JsonSlurper
import groovy.json.JsonParserType

/**
 * Check if a tool should run based on the tool_configs JSON string.
 *
 * @param json_text JSON string from PARSE_FRAGPIPE_WORKFLOW tool_configs channel
 * @param tool_name Name of the tool to check (e.g., 'msfragger', 'msbooster')
 * @return true if the tool should run, false otherwise
 */
def shouldRunTool(json_text, tool_name) {
    def configs = new JsonSlurper()
        .setType(JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.run ?: false
}

/**
 * Get the args string for a tool from the tool_configs JSON string.
 *
 * @param json_text JSON string from PARSE_FRAGPIPE_WORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Args string for the tool, or empty string if not found
 */
def getToolArgs(json_text, tool_name) {
    def configs = new JsonSlurper()
        .setType(JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.args ?: ''
}

/**
 * Get the modmasses string for a tool from the tool_configs JSON string.
 *
 * Used by IonQuant for --modlist (MBR modification mass tracking).
 * Returns comma-separated modification masses, or empty string if not set.
 *
 * @param json_text JSON string from PARSE_FRAGPIPE_WORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Comma-separated modmasses string, or empty string if not found
 */
def getToolModmasses(json_text, tool_name) {
    def configs = new JsonSlurper()
        .setType(JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.modmasses ?: ''
}

/**
 * Get an arbitrary field value for a tool from the tool_configs JSON string.
 *
 * Used to extract individual algorithm parameters that the parser emits as
 * separate JSON fields (e.g., FPOP's region_size, control_label).
 *
 * @param json_text JSON string from PARSE_FRAGPIPE_WORKFLOW tool_configs channel
 * @param tool_name Name of the tool (e.g., 'fpop')
 * @param field_name Name of the field to retrieve (e.g., 'region_size')
 * @param default_val Default value if the field is not found
 * @return Field value as a string, or default_val if not found
 */
def getToolField(json_text, tool_name, field_name, default_val = '') {
    def configs = new JsonSlurper()
        .setType(JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.get(field_name) ?: default_val
}

/**
 * Get the report_args string for a tool from the tool_configs JSON string.
 *
 * Used by philosopher filter to pass report flags (--decoys, --removecontam)
 * to the `philosopher report` command.
 *
 * @param json_text JSON string from PARSE_FRAGPIPE_WORKFLOW tool_configs channel
 * @param tool_name Name of the tool
 * @return Report args string, or empty string if not found
 */
def getToolReportArgs(json_text, tool_name) {
    def configs = new JsonSlurper()
        .setType(JsonParserType.CHARACTER_SOURCE)
        .parseText(json_text)
    return configs?.get(tool_name)?.report_args ?: ''
}

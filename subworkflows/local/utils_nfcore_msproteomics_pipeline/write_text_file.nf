/*
 * Helper of WRITETEXTFILE (modules/local/writetextfile), shared by refactors A and B:
 * the bash that writes a file body byte for byte. Unit tests in tests/write_text_file.nf.test.
 */

/**
 * Bash that writes a file body byte for byte with a quoted here-doc (WRITETEXTFILE).
 *
 * @param name File name
 * @param text File body; a body without a final newline is written without one
 * @return Bash command
 */
def textFileCommand(name, text) {
    def body = text.toString()
    def delimiter = 'WRITE_TEXT_FILE_END'
    def quoted = "'" + name.toString().replace("'", "'\\''") + "'"
    if (body.readLines().contains(delimiter)) {
        error("WRITETEXTFILE: ${name} contains the here-doc delimiter line ${delimiter}")
    }
    if (body.isEmpty()) {
        return ": > ${quoted}".toString()
    }
    if (body.endsWith('\n')) {
        return "cat > ${quoted} <<'${delimiter}'\n${body.substring(0, body.length() - 1)}\n${delimiter}".toString()
    }
    return "printf '%s' \"\$(cat <<'${delimiter}'\n${body}\n${delimiter}\n)\" > ${quoted}".toString()
}

/*
 * Test-only helpers for tests/test_data_base_path.nf.test.
 */

/**
 * Remote spectra URLs in the test samplesheets that do not start with the test data base path.
 *
 * @param samplesheet_dir Directory holding the test_*.csv samplesheets
 * @param base_path params.pipelines_testdata_base_path
 * @return List of "<file>: <url>" entries outside base_path (empty when consistent)
 */
def remoteUrlsOutsideBase(samplesheet_dir, base_path) {
    def offending = []
    new File(samplesheet_dir.toString()).listFiles().findAll { f -> f.name ==~ /test_.*\.csv/ }.sort { f -> f.name }.each { f ->
        f.readLines().drop(1).each { line ->
            line.split(',').findAll { cell -> cell ==~ /^[a-z0-9]+:\/\/.*/ && !cell.startsWith(base_path + '/') }.each { cell -> offending << "${f.name}: ${cell}" }
        }
    }
    return offending
}

/**
 * Test-data URLs in a test file that do not start with the test-data base path.
 *
 * @param test_file nf-test file to scan
 * @param base_path params.pipelines_testdata_base_path
 * @return "line: url" entries outside the base path
 */
def testDataUrlsOutsideBase(test_file, base_path) {
    def offending = []
    new File(test_file.toString()).readLines().eachWithIndex { line, i ->
        (line =~ /https:\/\/raw\.githubusercontent\.com\/[^'"\s]*test-datasets[^'"\s]*/).each { url ->
            if (!url.startsWith(base_path + '/')) {
                offending << "${i + 1}: ${url}"
            }
        }
    }
    return offending
}

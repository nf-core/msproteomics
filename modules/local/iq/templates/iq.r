#!/usr/bin/env Rscript

# if diann doesn't exist, install it
if (!requireNamespace("iq", quietly = TRUE)) {
    install.packages('iq', dependencies = TRUE, repos='http://cran.rstudio.com/')
}

library(iq)
library(readr)

# Load a DIA-NN report (a small sample report is included in this repository)

report <- file.path("${report_path}")
q <- as.numeric("${q}")
pg.q <- as.numeric("${pgq}")
contaminant_pattern <- "${contaminant_pattern}"

SAMPLE_ID <- "Run"
PRIMARY_ID <- "Protein.Group"
ANNOTATION_COLS <- c("Protein.Names", "Genes")
OUTPUT_FILENAME <- "maxlfq.tsv"

df <- read_tsv(report)
if (nrow(df) == 0) {
    stop(paste0("ERROR: IQ: DIA-NN report '", report, "' has no rows; DIA-NN identified nothing"))
}
if (contaminant_pattern != "" && !is.null(contaminant_pattern)) {
    # Remove contaminant proteins
    df <- df[!grepl(contaminant_pattern, df\$Genes),]
    report <- paste0("contaminants_removed_", report)
    write_tsv(df, report)
    if (nrow(df) == 0) {
        stop(paste0("ERROR: IQ: every row of the DIA-NN report matches contaminant_pattern '", contaminant_pattern, "'"))
    }
}

# The same thresholds as filter_double_less below; iq keeps a row when value < threshold
fdr_filter <- c(
    "Q.Value" = q,
    "PG.Q.Value" = pg.q,
    "Lib.Q.Value" = q,
    "Lib.PG.Q.Value" = pg.q)
passing <- Reduce(`&`, lapply(names(fdr_filter), function(col) df[[col]] < fdr_filter[[col]]))
filter_desc <- paste(paste0(names(fdr_filter), " < ", fdr_filter), collapse = " & ")


# Protein x samples MaxLFQ matrix filtered at precursor and protein-group FDR.
# When no row passes (e.g. DIA-NN reports PG.Q.Value = 1 for every row because a small
# search space gives too few protein groups to estimate protein-group FDR),
# process_long_format() reports "# quantitative values after filtering = 0", returns
# without writing the output file, and R exits 0. Write a header-only table instead.
if (!any(passing)) {
    message(paste0("WARNING: IQ: ", OUTPUT_FILENAME, ": no report row passes ", filter_desc,
                   " (", nrow(df), " rows; min Q.Value = ", min(df\$Q.Value),
                   ", min PG.Q.Value = ", min(df\$PG.Q.Value),
                   ", min Lib.Q.Value = ", min(df\$Lib.Q.Value),
                   ", min Lib.PG.Q.Value = ", min(df\$Lib.PG.Q.Value), "); writing a header-only table"))
    writeLines(paste(c(PRIMARY_ID, ANNOTATION_COLS, unique(df[[SAMPLE_ID]])), collapse = "\\t"), OUTPUT_FILENAME)
} else {
    process_long_format(
        report,
        sample_id = SAMPLE_ID,
        primary_id = PRIMARY_ID,
        intensity_col = "Fragment.Quant.Raw",
        output_filename = OUTPUT_FILENAME,
        annotation_col = ANNOTATION_COLS,
        filter_double_less = fdr_filter
    )
    if (!file.exists(OUTPUT_FILENAME)) {
        stop(paste0("ERROR: IQ: ", sum(passing), " report rows pass ", filter_desc,
                    " but process_long_format() wrote no ", OUTPUT_FILENAME))
    }
}

#
# save versions file
#
versions_file <- file("versions.yml")
write(
    paste(
        '${task.process}:',
        paste0('  r-base: "', R.Version()\$version.string, '"'),
        paste0('  iq: "', as.character(packageVersion("iq")), '"'),
        sep = "\\n"
    ),
    versions_file
)
close(versions_file)

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

df <- read_tsv(report)
if (contaminant_pattern != "" && !is.null(contaminant_pattern)) {
    # Remove contaminant proteins
    df <- df[!grepl(contaminant_pattern, df\$Genes),]
    report <- paste0("contaminants_removed_", report)
    write_tsv(df, report)
}


# Precursors x samples matrix filtered at 1% precursor and protein group FDR
process_long_format(
    report,
    sample_id = "Run",
    intensity_col = "Fragment.Quant.Raw",
    output_filename = "maxlfq.tsv",
    annotation_col = c("Protein.Names", "Genes"),
    filter_double_less = c(
        "Q.Value" = q,
        "PG.Q.Value" = pg.q,
        "Lib.Q.Value" = q,
        "Lib.PG.Q.Value" = pg.q)
)

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

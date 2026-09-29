#!/usr/bin/env Rscript

# if diann doesn't exist, install it
if (!requireNamespace("diann", quietly = TRUE)) {
    devtools::install_github('https://github.com/vdemichev/diann-rpackage')
}

library(diann)
library(readr)

# Load a DIA-NN report (a small sample report is included in this repository)

report <- file.path("${report_path}")
q <- as.numeric("${q}")
pg.q <- as.numeric("${pgq}")
contaminant_pattern <- "${contaminant_pattern}"

# diann_matrix() filters at 1% precursor FDR unless q is given (see the calls below)
DIANN_MATRIX_DEFAULT_Q <- 0.01
SAMPLE_HEADER <- "File.Name"

df <- diann_load(report)
if (nrow(df) == 0) {
    stop(paste0("ERROR: DIANNR: DIA-NN report '", report, "' has no rows; DIA-NN identified nothing"))
}
if (contaminant_pattern != "" && !is.null(contaminant_pattern)) {
    # Remove contaminant proteins
    df <- df[!grepl(contaminant_pattern, df\$Genes),]
    report <- paste0("contaminants_removed_", report)
    write_tsv(df, report)
    if (nrow(df) == 0) {
        stop(paste0("ERROR: DIANNR: every row of the DIA-NN report matches contaminant_pattern '", contaminant_pattern, "'"))
    }
}
samples <- unique(df[[SAMPLE_HEADER]])

# Write a quantity matrix. When no report row passes the FDR filter (e.g. DIA-NN
# reports PG.Q.Value = 1 for every row because a small search space gives too few
# protein groups to estimate protein-group FDR), write a header-only table (the sample
# columns) and say so. An empty result from rows that did pass is a bug: fail.
write_quant_table <- function(quant, path, n_passing, filter_desc) {
    if (n_passing == 0) {
        message(paste0("WARNING: DIANNR: ", path, ": no report row passes ", filter_desc,
                       " (", nrow(df), " rows; min Q.Value = ", min(df\$Q.Value),
                       ", min PG.Q.Value = ", min(df\$PG.Q.Value), "); writing a header-only table"))
        writeLines(paste(samples, collapse = "\t"), path)
        return(invisible(NULL))
    }
    if (is.null(quant) || nrow(quant) == 0) {
        stop(paste0("ERROR: DIANNR: ", path, ": ", n_passing, " report rows pass ", filter_desc,
                    " but the quantity table is empty"))
    }
    write.table(quant, file.path(path), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}

pass_matrix <- df\$Q.Value <= DIANN_MATRIX_DEFAULT_Q & df\$PG.Q.Value <= pg.q
matrix_filter_desc <- paste0("Q.Value <= ", DIANN_MATRIX_DEFAULT_Q, " & PG.Q.Value <= ", pg.q)

# Precursors x samples matrix filtered at 1% precursor and protein group FDR
precursors <- diann_matrix(df, pg.q = pg.q)
write_quant_table(precursors, "precursors.tsv", sum(pass_matrix), matrix_filter_desc)

# Peptides without modifications - taking the maximum of the respective precursor quantities
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = pg.q)
write_quant_table(peptides, "peptides.tsv", sum(pass_matrix), matrix_filter_desc)

# Peptides without modifications - using the MaxLFQ algorithm.
# diann_maxlfq() cannot take an empty table: its loop runs over 1:length(groups) = 1:0
# and fails with "subscript out of bounds", so it is only called when rows pass.
pass_peptide_lfq <- df\$Q.Value <= q
peptides.maxlfq <- NULL
if (any(pass_peptide_lfq)) {
    peptides.maxlfq <- diann_maxlfq(df[pass_peptide_lfq,], id.header = "Stripped.Sequence", quantity.header = "Precursor.Normalised")
}
write_quant_table(peptides.maxlfq, "peptides_maxlfq.tsv", sum(pass_peptide_lfq), paste0("Q.Value <= ", q))

# Genes identified and quantified using proteotypic peptides
unique.genes <- diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = pg.q)
write_quant_table(unique.genes, "unique_genes.tsv", sum(pass_matrix & df\$Proteotypic != 0),
                  paste0(matrix_filter_desc, " & Proteotypic"))

# Protein group quantities using MaxLFQ algorithm
pass_pg_lfq <- df\$Q.Value <= pg.q & df\$PG.Q.Value <= pg.q
protein.groups <- NULL
if (any(pass_pg_lfq)) {
    protein.groups <- diann_maxlfq(df[pass_pg_lfq,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
}
write_quant_table(protein.groups, "protein_groups_maxlfq.tsv", sum(pass_pg_lfq),
                  paste0("Q.Value <= ", pg.q, " & PG.Q.Value <= ", pg.q))

#
# save versions file
#
versions_file <- file("versions.yml")
write(
    paste(
        '${task.process}:',
        paste0('  r-base: "', R.Version()\$version.string, '"'),
        paste0('  diann: "', as.character(packageVersion("diann")), '"'),
        sep = "\\n"
    ),
    versions_file
)
close(versions_file)

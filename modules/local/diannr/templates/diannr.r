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

df <- diann_load(report)
if (contaminant_pattern != "" && !is.null(contaminant_pattern)) {
    # Remove contaminant proteins
    df <- df[!grepl(contaminant_pattern, df\$Genes),]
    report <- paste0("contaminants_removed_", report)
    write_tsv(df, report)
}


# Precursors x samples matrix filtered at 1% precursor and protein group FDR
precursors <- diann_matrix(df, pg.q = pg.q)

write.table(precursors , file.path("precursors.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Peptides without modifications - taking the maximum of the respective precursor quantities
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = pg.q)

write.table(peptides , file.path("peptides.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Peptides without modifications - using the MaxLFQ algorithm
peptides.maxlfq <- diann_maxlfq(df[df\$Q.Value <= q,], id.header = "Stripped.Sequence", quantity.header = "Precursor.Normalised")
write.table(peptides.maxlfq , file.path("peptides_maxlfq.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Genes identified and quantified using proteotypic peptides
unique.genes <- diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = pg.q)

write.table(unique.genes , file.path("unique_genes.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


# Protein group quantities using MaxLFQ algorithm
protein.groups <- diann_maxlfq(df[df\$Q.Value <= pg.q & df\$PG.Q.Value <= pg.q,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

write.table(protein.groups , file.path("protein_groups_maxlfq.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

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

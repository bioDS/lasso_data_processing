#!/usr/bin/env Rscript

require(biomaRt)
require(stringr)
require(data.table)
require(normInfectX)

test_sirna = vaccinia[25435,]

# some tests
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

test = getBM(attributes = "entrezgene_id", filters="refseq_mrna", values=c("nm_021222", "nm_001350721"), mart=ensembl)

sequence = as.character(test_sirna$Sequence_antisense_5_3)
name = as.character(test_sirna$Catalog_number)
fasta_file = sprintf("%s.fasta", name)
output_file = sprintf("risearch_%s.out.gz", name)

write(sprintf("> %s", name), file=fasta_file)
write(sequence, file=fasta_file, append=TRUE)

system(sprintf("./RIsearch-2.1/bin/risearch2.x --index=mrna.suf --query=%s --seed=2:8 --energy=-25", fasta_file))

dt = fread(output_file)
dt$mrna = str_extract(dt$V4, "[^.]*")

suppressed_genes = getBM(attributes = "entrezgene_id", filters="refseq_mrna", values=dt$mrna, mart=ensembl)


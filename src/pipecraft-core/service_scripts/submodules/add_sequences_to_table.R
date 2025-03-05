#!/usr/bin/env Rscript

# Prepare OTU tables for merging sequencing runs
 # 1. add sequences to the table
 # 2. export sequences with the size annotation

# load variables
args <- commandArgs(trailingOnly = TRUE)
otu_table_file <- args[1]
fasta_file <- args[2]

### Load OTU table
suppressMessages(library(data.table))
cat(";; Loading input table:", otu_table_file, "\n")
otu_table = fread(
  file = otu_table_file,
  sep = "\t", header = TRUE)

base_name = tools::file_path_sans_ext(basename(otu_table_file))
output_name = paste0(base_name, "_wSeqs.txt")

### Load fasta file and add sequences back to the table
cat(";; Adding sequences back to the table\n")
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr)) 
# read the FASTA file
cat(";; Loading input FASTA:", fasta_file, "\n")
fasta_sequences = readDNAStringSet(fasta_file)
sequences = as.character(fasta_sequences)
names(sequences) = names(fasta_sequences)

# add the sequences as 2nd column
otu_table = otu_table %>%
  mutate(Sequence = sequences[OTU]) %>%
  select(OTU, Sequence, everything())

## Export table
cat(";; Exporting table with sequences (2nd col is sequence) \n")
cat(";; ", output_name, "\n")
fwrite(x = otu_table,
  file = output_name,
  sep = "\t", compress = "none")

cat(";; Done \n")

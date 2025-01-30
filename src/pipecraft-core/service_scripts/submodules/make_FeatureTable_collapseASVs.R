#!/usr/bin/Rscript

## Script to summarize ASV abundance by sample (within each collapsed ASV based on usearch_global)
# Used to collapse ASVs in the DADA2 workflow (replacing "collapseNoMismatch" function of dada2)
# Edited 22.01.2025, form the ASV_OTU_merging_script.R
  # 30.01.2025: added sequences to the output table (for merge_runs); reads fasta file

## Usage:
# ./make_FeatureTable_collapseASVs.R \
#    --derepuc  "ASVs_derep.uc" \
#    --uc       "usearch_global.uc" \
#    --asv      "ASV_table_long.txt" \
#    --fasta    "sequences.fasta" \
#    --output   "ASVs_table_collapsed.txt"

## Input:
# ASVs_derep.uc            = UC file from `vsearch --derep_fulllenth ...`
# usearch_global.uc        = UC file from `vsearch --usearch_global ... -uc usearch_global.uc`
# ASV_table_long.txt       = tab-delimited 3-column table (SequenceID, SampleID, Abundance), with header
# sequences.fasta           = FASTA file with sequences
# ASVs_table_collapsed.txt = resulting file

suppressMessages(library(optparse))
suppressMessages(library(data.table))

## Parse arguments
option_list <- list(
  make_option(c("-d", "--derepuc"), action="store", default=NA, type='character', help="UC-file from the dereplication step"),
  make_option(c("-u", "--uc"),      action="store", default=NA, type='character', help="UC-file from the clustering step"),
  make_option(c("-a", "--asv"),     action="store", default=NA, type='character', help="ASV sequence counts per sample"),
  make_option(c("-f", "--fasta"),   action="store", default=NA, type='character', help="FASTA file, the one used for the dereplication(collapse) step"),
  make_option(c("-o", "--output"),  action="store", default="ASVs_table_collapsed.txt", type='character', help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

INP_DEREPUC <- opt$derepuc       # e.g., INP_DEREPUC <- "ASVs_derep.uc"
INP_UC      <- opt$uc            # e.g., INP_UC  <- "usearch_global.uc"
INP_ASV     <- opt$asv           # e.g., INP_ASV <- "ASV_table_long.txt"
INP_FASTA   <- opt$fasta         # e.g., INP_FASTA <- "ASVs.fasta (the one used for the dereplication(collapse) step)"
OUTPUT      <- opt$output        # e.g., OUTPUT <- "ASVs_table_collapsed.txt"

## Load input data - ASV table
cat(";; Loading ASV table. File = ", INP_ASV, "\n")
ASV <- fread(file = INP_ASV, header = T, sep = "\t")
colnames(ASV) <- c("SequenceID", "SampleID", "Abundance")

## Load dereplication UC file
cat(";; Loading dereplication UC file. File = ", INP_DEREPUC, "\n")
DUC <- fread(file = INP_DEREPUC, header = FALSE, sep = "\t")
DUC <- DUC[ V1 != "S" ]
DUC[, ASVall  := tstrsplit(V9, ";", keep = 1) ]
DUC[, ASVuniq := tstrsplit(V10, ";", keep = 1) ]
DUC[V1 == "C", ASVuniq := ASVall ]
DUC <- DUC[, .(ASVall, ASVuniq)]
DUC <- unique(DUC)

## Add dereplicated-ASV-IDs to the ASV table
cat(";; Adding dereplicated-ASV-IDs to the ASV table.\n")
ASV <- merge(x = ASV, y = DUC,
  by.x = "SequenceID", by.y = "ASVall",
  all.x = TRUE)

## Clean up, to save some RAM
rm(DUC)
gc()

## Load clustering UC file
cat(";; Loading clustering UC file. File = ", INP_UC, "\n")
UC <- fread(file = INP_UC, header = FALSE, sep = "\t")
UC <- UC[ V1 != "S" ]
UC[, ASV := tstrsplit(V9, ";", keep = 1) ]
UC[, OTU := tstrsplit(V10, ";", keep = 1) ]
UC[V1 == "C", OTU := ASV ]
UC <- UC[, .(ASV, OTU)]

## Add OTU IDs to the ASV table
cat(";; Adding OTU IDs to the ASV table.\n")
ASV <- merge(x = ASV, y = UC,
  by.x = "ASVuniq", by.y = "ASV",
  all.x = TRUE)

## Clean up, to save some RAM
rm(UC)
gc()
ASV[, ASVuniq := NULL ]

## Check if there are some missing data
cat(";; Checking if there are some missing data.\n")
if(any(is.na(ASV$OTU))){
  nreads   <- sum(ASV[ is.na(OTU) ]$Abundance)
  nrecords <- sum(is.na(ASV$OTU))
  cat("WARNING: ", nreads, " reads (in ", 
    nrecords, "records) were not assigned to OTUs\n")
  ASV <- ASV[ !is.na(OTU) ]
}

## Summarize ASV abundance
cat(";; Summarizing ASV abundance.\n")
OTU <- ASV[ , 
  .(Abundance = sum(Abundance, na.rm = TRUE)), 
  by = c("SampleID", "OTU")]

## Reshape to wide format
RES <- dcast(data = ASV,
  formula = OTU ~ SampleID,
  value.var = "Abundance",
  fun.aggregate = sum, fill = 0)

## Sort by total OTU/ASV abundance
cat(";; Sorting by total OTU/ASV abundance.\n")
TotAb <- rowSums(x = RES[, -1])
RES <- RES[ order(TotAb, decreasing = T) , ]

## Add sequences back to the table
cat(";; Adding sequences back to the table\n")
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
# read the FASTA file
cat(";; Loading input FASTA:", INP_FASTA, "\n")
fasta_sequences = readDNAStringSet(INP_FASTA)
sequences = as.character(fasta_sequences)
names(sequences) = names(fasta_sequences)
# add the sequences as 2nd column
RES = RES %>%
  mutate(Sequence = sequences[OTU]) %>%
  select(OTU, Sequence, everything())

## Export table
cat(";; Exporting table (2nd col is sequence)\n")
cat(OUTPUT, "\n")
fwrite(x = RES,
  file = OUTPUT,
  sep = "\t", compress = "none")

cat(";; Done\n")

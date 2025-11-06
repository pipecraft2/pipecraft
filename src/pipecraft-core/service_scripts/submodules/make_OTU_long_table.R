#!/usr/bin/Rscript

## Make OTU table based on dereplicated sequences and UC file
 # This is the same as ASV_OTU_merging_script.R, but outputs OTU table in long format

## Usage:
# ./make_OTU_table_multiRunDir.R \
#    --derepuc  "Glob_derep.uc" \
#    --uc       "OTUs.uc" \
#    --asv      "ASV_table_long.txt" \
#    --rmsingletons TRUE \
#    --output   "OTU_table_long.txt"

## Input:
# Glob_derep.uc = UC file from `vsearch --derep_fulllenth ...`
# OTUs.uc       = UC file from `vsearch --cluster_size ... -uc OTUs.uc`
# ASV_table_long.txt = tab-delimited 3-column table (SampleID, SequenceID, Abundance), no header
# OTU_table_long.txt = resulting file

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

## Parse arguments
option_list <- list(
  make_option(c("-d", "--derepuc"), action="store", default=NA, type='character', help="UC-file from the dereplication step"),
  make_option(c("-u", "--uc"),      action="store", default=NA, type='character', help="UC-file from the clustering step"),
  make_option(c("-a", "--asv"),     action="store", default=NA, type='character', help="ASV sequence counts per sample"),
  make_option(c("-s", "--rmsingletons"), action="store", default=TRUE, type='logical', help="Remove global singletons"),
  make_option(c("-f", "--fasta"),   action="store", default=NA, type='character', help="FASTA file containing OTU sequences"),
  make_option(c("-o", "--output"),  action="store", default="OTU_table.txt", type='character', help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

INP_DEREPUC <- opt$derepuc       # e.g., INP_DEREPUC <- "Glob_derep.uc"
INP_UC      <- opt$uc            # e.g., INP_UC  <- "OTUs.uc"
INP_ASV     <- opt$asv           # e.g., INP_ASV <- "ASV_table_long.txt"
INP_FASTA   <- opt$fasta         # e.g., INP_FASTA <- "OTUs.fasta"
RMSINGLETON <- opt$rmsingletons  # e.g., RMSINGLETON <- TRUE
OUTPUT      <- opt$output        # e.g., OUTPUT <- "OTU_table.txt"

cat(";; Making OTU table ...\n")
## Load input data - ASV table
ASV <- fread(file = INP_ASV, header = FALSE, sep = "\t",
  col.names = c("SampleID", "SequenceID", "Abundance"))

## Load dereplication UC file
DUC <- fread(file = INP_DEREPUC, header = FALSE, sep = "\t")
DUC <- DUC[ V1 != "S" ]
DUC[, ASVall  := tstrsplit(V9, ";", keep = 1) ]
DUC[, ASVuniq := tstrsplit(V10, ";", keep = 1) ]
DUC[V1 == "C", ASVuniq := ASVall ]
DUC <- DUC[, .(ASVall, ASVuniq)]
DUC <- unique(DUC)

## Add dereplicated-ASV-IDs to the ASV table
ASV <- merge(x = ASV, y = DUC,
  by.x = "SequenceID", by.y = "ASVall",
  all.x = TRUE)

## Clean up, to save some RAM
rm(DUC)
invisible(gc())

## Load clustering UC file
UC <- fread(file = INP_UC, header = FALSE, sep = "\t")
UC <- UC[ V1 != "S" ]
UC[, ASV := tstrsplit(V9, ";", keep = 1) ]
UC[, OTU := tstrsplit(V10, ";", keep = 1) ]
UC[V1 == "C", OTU := ASV ]
UC <- UC[, .(ASV, OTU)]

## Add OTU IDs to the ASV table
ASV <- merge(x = ASV, y = UC,
  by.x = "ASVuniq", by.y = "ASV",
  all.x = TRUE)

## Clean up, to save some RAM
rm(UC)
ASV[, ASVuniq := NULL ]

## Check if there are some missing data
if(any(is.na(ASV$OTU))){
  nreads   <- sum(ASV[ is.na(OTU) ]$Abundance)
  nrecords <- sum(is.na(ASV$OTU))
  cat("WARNING: ", nreads, " reads (in ", 
    nrecords, "records) were not assigned to OTUs\n")
  ASV <- ASV[ !is.na(OTU) ]
}

## Remove global singletons
if(RMSINGLETON){

  ## Estimate total OTU abundance
  OTUab <- ASV[, .(TotalAbundance = sum(Abundance)), by = "OTU" ]
  SINGLETONS <- OTUab[ TotalAbundance <= 1 ]
  
  if(nrow(SINGLETONS) > 0){
    cat("Removing ", nrow(SINGLETONS), " singleton OTUs\n")
    ASV <- ASV[ ! OTU %in% SINGLETONS$OTU ]
  }
}

## Summarize ASV abundance
OTU <- ASV[ , 
  .(Abundance = sum(Abundance, na.rm = TRUE)), 
  by = c("SampleID", "OTU")]

## Sort by OTU name
RES <- OTU[ order(OTU, SampleID) , ]

invisible(gc())

## Export table in long format
cat(";; Exporting table in long format\n")
cat(OUTPUT, "\n")
fwrite(x = RES,
  file = OUTPUT,
  sep = "\t", compress = "none")

cat(";; Done\n")

#!/usr/bin/env Rscript

## Script to perform tag-jump removal

## To do:
#  - add HMM profile ID if ITSx was used

# Input is given as positional arguments:
#   1. non-filtered Seq table   (`Seq_tab_not_filtered.txt.gz`)
#   2. Sequences in fasta       (`Seq_not_filtered.fa.gz`)
#   3. sequence mapping to OTUs (`Sample_mapping.uc.gz`)
#   4. tag-jumped OTU list      (`TagJump_OTUs.RData`)
#   5. de novo chimera scores   (`DeNovo_Chimera.txt`)
#   6. sequence qualities       (`SeqQualities.parquet`)

# Outputs:
#  - FASTA with filtered Seqs       `Seqs.fa.gz`
#  - Seq table in long format       `Seqs.txt.gz`  (with additional sequence info)
#  - Data in Parquet format         `Seqs.parquet`
#  - Seq table in wide format       `Seq_tab.txt.gz` -- deprecated
#  - Data in R-serialization format `Seqs.RData`     -- deprecated


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")

load_pckg("optparse")
load_pckg("data.table")
load_pckg("Biostrings")
load_pckg("plyr")
load_pckg("arrow")
# load_pckg("dplyr")
# load_pckg("openxlsx")


cat("Parsing input options and arguments...\n")

option_list <- list(
  make_option("--seqtab",  action="store", default=NA, type='character', help = "Non-filtered sequence table (tab-delimited)"),
  make_option("--fasta",   action="store", default=NA, type='character', help = "Sequences in FASTA format"),
  make_option("--mapping", action="store", default=NA, type='character', help = "Sequence mapping to OTUs (UC format)"),
  make_option("--tagjump", action="store", default=NA, type='character', help = "Tag-jumped OTU list (RData format)"),
  make_option("--chimera", action="store", default=NA, type='character', help = "De novo chimera scores"),
  make_option("--quality", action="store", default=NA, type='character', help = "Sequence qualities (Parquet format)"),
  make_option("--threads", action="store", default=4,  type='integer',   help = "Number of CPU threads to use")
)

opt <- parse_args(OptionParser(option_list=option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Validation of the required arguments
required_args <- c("seqtab", "fasta", "mapping", "tagjump", "quality")
missing_args <- required_args[ sapply(required_args, function(x) is.na(opt[[x]])) ]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse=", "))
}

## Assign variables
SEQTAB  <- opt$seqtab
FASTA   <- opt$fasta
MAPPING <- opt$mapping
TAGJUMP <- opt$tagjump
CHIMERA <- opt$chimera
QUALITY <- opt$quality
CPUTHREADS <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("Non-filtered sequence table: ",  SEQTAB,     "\n", sep=""))
cat(paste("Sequences in FASTA format: ",    FASTA,      "\n", sep=""))
cat(paste("Sequence mapping to OTUs: ",     MAPPING,    "\n", sep=""))
cat(paste("Tag-jumped OTU list: ",          TAGJUMP,    "\n", sep=""))
cat(paste("De novo chimera scores: ",       CHIMERA,    "\n", sep=""))
cat(paste("Sequence qualities: ",           QUALITY,    "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS, "\n", sep=""))

cat("\n")


## Debug:
# SEQTAB     <- "Seq_tab_not_filtered.txt.gz"
# FASTA      <- "Seq_not_filtered.fa.gz"
# MAPPING    <- "Sample_mapping.uc.gz"
# TAGJUMP    <- "TagJump_OTUs.RData"
# CHIMERA    <- "DeNovo_Chimera.txt"
# QUALITY    <- "SeqQualities.parquet"
# CPUTHREADS <- 4


## Set CPU thread number
cat("\nSetting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)     # for data.table
set_cpu_count(CPUTHREADS)              # for arrow

######################################
###################################### Load the data
######################################

## Load sequnece table
cat("\n\n..Loading non-filtered sequence table\n")

TAB <- fread(
  file = SEQTAB,
  sep = "\t", header = FALSE,
  col.names = c("SeqID", "Abundance", "SampleID"))


## Load sequences in fasta format
cat("..Loading sequneces in FASTA format\n")
SQS <- readDNAStringSet(filepath = FASTA)


## Load sequence mapping to pre-clustered groups for tag-jump removal
cat("..Loading sequence mapping table\n")
MAP <- fread(
  file = MAPPING,
  header = FALSE, sep = "\t")

MAP <- MAP[ V1 != "S" ]
MAP[, c("SeqID", "SampleID") := tstrsplit(V9, ";", keep = c(1,3)) ]
MAP[, OTU := tstrsplit(V10, ";", keep = 1) ]
MAP[V1 == "C", OTU := SeqID ]
MAP[, SampleID := gsub(pattern = "sample=", replacement = "", x = SampleID) ]
MAP <- MAP[, .(SeqID, SampleID, OTU) ]


## Load list of tag-jumped OTUs
cat("..Loading list of tag-jumped OTUs\n")
JMP <- readRDS( TAGJUMP )


## Load de novo chimera scores
cat("..Loading de novo chimera scores\n")

CHI <- try(
  fread(
    file = CHIMERA,
    header = FALSE, sep = "\t",
    col.names = c("SeqID", "DeNovo_Chimera_Score", "SampleID"))
  )

if("try-error" %in% class(CHI)){
  cat("\nCould not read the file with de novo chimeric scores\n")
  cat("Most likely, the file file is empty (no de novo chimeras)\n")
  
  ## Initialize empty data table
  CHI <- data.table(SeqID = character(), DeNovo_Chimera_Score = numeric(), SampleID = character())
}


## Load sequence quality scores
cat("..Loading sequence quality scores\n")
QLT <- arrow::open_dataset(QUALITY) |>
  dplyr::select(Hash, Length, AvgPhredScore, MaxEE, MEEP) |>
  dplyr::collect() |>
  dplyr::filter(Hash %in% unique(TAB$SeqID)) |>
  setDT()

setnames(QLT,
  old = c("Hash", "Length", "AvgPhredScore"),
  new = c("SeqID", "SeqLen", "PhredScore"))

## Quality data:
# old header: c("SampleID", "SeqID", "SeqLen", "PhredScore", "MaxEE", "MEEP")
# new header: c("SampleID", "Hash", "PacBioID", "PhredScore", "MaxEE", "MEEP", "Sequence", "Quality", "Length")


## Create SeqID___SampleID column
TAB[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
# QLT[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
MAP[, SeqID___SampleID := paste0(SeqID, "___", SampleID) ]

MAP[, c("SeqID", "SampleID") := NULL ]

######################################
###################################### Remove tag-jumps
######################################

cat("\n\n..Removing tag-jumped sequences\n")

if(nrow(JMP) > 0){

  # JMP[ , SeqID___SampleID := paste0(OTU, "___", SampleID) ]
  JMP[ , TagJump := TRUE ]

  ## Add OTU ID to sequences
  cat("...Adding OTU IDs to sequences\n")

  TAB <- merge(x = TAB, y = MAP,
    by = "SeqID___SampleID", all.x = TRUE)

  ## Add tag-jump information to the sequence table
  TAB <- merge(x = TAB, y = JMP,
    by = c("OTU", "SampleID"),
    all.x = TRUE)

  cat("... ", sum(TAB$TagJump, na.rm = TRUE), " tag-jump occurrences detected\n")

  ## Remove tag-jumps
  if(any(TAB$TagJump)){
    TAB <- TAB[(!TagJump %in% TRUE)]
    # because of NAs, TAB[ ! TagJump ] does not work properly
  }

  TAB[, TagJump := NULL ]
  TAB[, OTU := NULL ]

# end of `nrow(JMP) > 0`
} else {

  cat("...no tag-jumps found\n")

}


######################################
###################################### Add quality scores, for non-singleton use max score
######################################

cat("\n\n..Adding quality scores\n")

cat("...Prepareing quality scores\n")
setorder(QLT, SeqID, -PhredScore)
QLT <- QLT[QLT[, .I[which.max(PhredScore)], by=SeqID]$V1]

cat("...Adding data to the main table\n")
if(any(! TAB$SeqID %in% QLT$SeqID)){
  cat("WARNING: Some sequences are not present in the quality table\n")
}

TAB <- merge(x = TAB, y = QLT, by = "SeqID", all.x = TRUE)

# with(TAB, plot(Abundance, PhredScore))


######################################
###################################### Add chimera info
######################################

cat("..Adding info about de novo chimeric sequences\n")

if(nrow(CHI) > 0){

  TAB <- merge(x = TAB, y = CHI,
    by = c("SeqID", "SampleID"), all.x = TRUE)

  ## Convert variables to numeric scores
  TAB[ , DeNovo_Chimera_Score := as.numeric(DeNovo_Chimera_Score) ]

  ## Classify sequences into putative chimeras
  TAB[ !is.na(DeNovo_Chimera_Score), DeNovo_Chimera := TRUE  ]
  TAB[  is.na(DeNovo_Chimera_Score), DeNovo_Chimera := FALSE ]

  cat("... ", sum( TAB$DeNovo_Chimera), " putative de novo chimeras found\n")
  cat("... ", sum(!TAB$DeNovo_Chimera), " non-chimeric sequences\n")

} else {

  ## No de novo chimeras

  TAB[ , DeNovo_Chimera_Score := as.numeric(NA) ]
  TAB[ , DeNovo_Chimera       := FALSE  ]

  cat("... ", 0,         " putative de novo chimeras found\n")
  cat("... ", nrow(TAB), " non-chimeric sequences\n")

}


######################################
###################################### Add sequences
######################################

cat("\n\n..Processing sequences\n")

SQTAB <- data.table(
  SeqHeader = names(SQS),
  Sequence = as.character(SQS))

## Split the header  (`feb76b9;size=1;sample=ABCD;` )
SQTAB[ , c("SeqID", "SampleID") := tstrsplit(x = SeqHeader, split = ";", keep = c(1,3)) ]
SQTAB[ , SeqHeader := NULL ]
SQTAB[ , SampleID := gsub(pattern = "sample=", replacement = "", x = SampleID) ]


SQTAB[ , SeqID___SampleID := paste0(SeqID, "___", SampleID) ]
SQTAB[ , c("SeqID", "SampleID") := NULL ]

cat("..Adding sequences to the main table\n")

TAB <- merge(x = TAB, y = SQTAB,
  by = c("SeqID___SampleID"), all.x = TRUE)


cat("..Sorting table by abundance, quality score, and SampleID\n")

setorder(x = TAB, -Abundance, -PhredScore, SampleID)


cat("..Preparing FASTA file with filtered sequences\n")

SQF <- DNAStringSet(x = TAB$Sequence)
names(SQF) <- paste0(TAB$SeqID, ";size=", TAB$Abundance, ";sample=", TAB$SampleID, ";")

## Export FASTA
cat("..Exporting FASTA file with filtered sequences\n")

writeXStringSet(x = SQF,
  filepath = "Seqs.fa.gz",
  compress = TRUE, format = "fasta", width = 9999)




######################################
###################################### Export results
######################################

# cat("..Reshaping sequence table into wide format\n")
# 
# TABW <- dcast(data = TAB,
#   formula = SeqID ~ SampleID, 
#   value.var = "Abundance",
#   fill = 0)


cat("..Exporting result\n")

setcolorder(
  x = TAB,
  neworder = c(
    "SeqID___SampleID", "SampleID", "SeqID",
    "Abundance", "SeqLen", "PhredScore", "MaxEE", "MEEP",
    "DeNovo_Chimera", "DeNovo_Chimera_Score",
    "Sequence"))

# cat("...Exporting RData\n")
# saveRDS(object = TAB, file = "Seqs.RData", compress = "xz")

cat("...Exporting Parquet\n")

write_parquet(
  x    = TAB,
  sink = "Seqs.parquet",
  compression       = "zstd",
  compression_level = 10,
  use_dictionary    = TRUE)


## Long table
cat("...Exporting long table\n")

TAB[ , SeqID___SampleID := NULL ]

fwrite(x = TAB, file = "Seqs.txt.gz", sep = "\t", compress = "gzip")

## Wide table
# cat("...Exporting wide table\n")

# fwrite(x = TABW, file = "Seq_tab.txt.gz", sep = "\t", compress = "gzip")


cat("All done.")

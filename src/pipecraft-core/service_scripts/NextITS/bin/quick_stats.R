#!/usr/bin/env Rscript

## Summarise number of reads (demultiplexed and primer-checked)

# quick_stats.R \
#   --raw          Counts_1.RawData.txt \
#   --qc           Counts_2.QC.txt \
#   --demuxed      Counts_3.Demux.txt \
#   --primer       Counts_4.PrimerCheck.txt \
#   --primerartef  Counts_4.PrimerArtefacts.txt \
#   --threads      4



############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(

  make_option("--raw",        action="store", default=NA, type='character', help="Raw read counts"),
  make_option("--qc",         action="store", default=NA, type='character', help="Counts of reads passed QC"),
  make_option("--demuxed",    action="store", default=NA, type='character', help="Counts of demultiplexed reads"),
  make_option("--primer",     action="store", default=NA, type='character', help="Counts of reads with both primers detected"),
  make_option("--primerartef",action="store", default=NA, type='character', help="Counts of primer artefacts"),
  make_option(c("-t", "--threads"), action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
# if(is.na(opt$raw)){
#   cat("Input file is not specified: ....\n", file=stderr())
#   stop()
# }



## Function to convert text "NA"s to NA
to_na <- function(x){ 
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Assign variables
RAW         <- opt$raw
QC          <- opt$qc
DEMUXED     <- opt$demuxed
PRIMER      <- opt$primer
PRIMERARTEF <- opt$primerartef
CPUTHREADS  <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("Counts - RawData: " ,     RAW, "\n", sep=""))
cat(paste("Counts - QC: " ,          QC, "\n", sep=""))
cat(paste("Counts - Demux: " ,       DEMUXED, "\n", sep=""))
cat(paste("Counts - PrimerCheck: " , PRIMER, "\n", sep=""))
cat(paste("Counts - Primer Artefacts: " ,   PRIMERARTEF, "\n", sep=""))
cat(paste("Number of CPU threads to use: ", CPUTHREADS,     "\n", sep=""))

cat("\n")


############################################## data for debuging

# RAW         <- "Counts_1.RawData.txt"
# QC          <- "Counts_2.QC.txt"
# DEMUXED     <- "Counts_3.Demux.txt"
# PRIMER      <- "Counts_4.PrimerCheck.txt"
# PRIMERARTEF <- "Counts_4.PrimerArtefacts.txt"
# CPUTHREADS  <- 6


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("metagMisc")
load_pckg("openxlsx")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table


######################################
###################################### Load the data
######################################

cat("\nLoading input data\n")


#### Per-dataset stats

## Load ASV table
cat("..Loading raw counts\n")
RAW <- fread(RAW)

cat("..Loading QC counts\n")
QC <- fread(QC)

#### Per-sample stats

SEQKITCOUNTS <- list()

cat("..Loading demux counts\n")
SEQKITCOUNTS$DEMUXED <- fread(DEMUXED)

cat("..Loading primer-checked data counts\n")
SEQKITCOUNTS$PRIMER      <- fread(PRIMER)
SEQKITCOUNTS$PRIMERARTEF <- fread(PRIMERARTEF)


## Remove NULL-files
null_seqk <- laply(.data = SEQKITCOUNTS, .fun = nrow)

if(any(null_seqk == 0)){
  cat("Some files with counts are missing:\n")
  to_rm <- which(null_seqk == 0)
  cat(".. ", paste(names(SEQKITCOUNTS)[ to_rm ], collapse = ", "), "\n")
  SEQKITCOUNTS[ to_rm ] <- NULL
  rm(to_rm)
}


## Process seqkit counts
seqkit_process <- function(x){
  if(nrow(x) > 0){

    ## Remove reudndant columns
    x <- x[ , .(file, num_seqs) ]

    ## Remove file extensions
    x[ , file := sub(pattern = ".fastq.gz",        replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fq.gz",           replacement = "", x = file) ]
    x[ , file := sub(pattern = ".fa.gz",           replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerChecked$",  replacement = "", x = file) ]
    x[ , file := sub(pattern = "_PrimerArtefacts$",replacement = "", x = file) ]
    
  }
  return(x)
}


cat("Processing data\n")
SEQKITCOUNTS <- llply(.data = SEQKITCOUNTS, .fun = seqkit_process)

## Rename columns
if(!is.null(SEQKITCOUNTS$DEMUXED)){
setnames(x = SEQKITCOUNTS$DEMUXED, old = "num_seqs", new = "Demultiplexed_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMER)){
setnames(x = SEQKITCOUNTS$PRIMER,  old = "num_seqs", new = "PrimerChecked_Reads", skip_absent = TRUE)
}
if(!is.null(SEQKITCOUNTS$PRIMERARTEF)){
setnames(x = SEQKITCOUNTS$PRIMERARTEF, old = "num_seqs", new = "PrimerArtefacts_Reads", skip_absent = TRUE)
}

## Merge seqkit and custom counts into a single list
cat("Pooling per-sample counts\n")
COUNTS <- SEQKITCOUNTS

## Pool per-file estimates
merge_dt <- function(x,y){ merge(x, y, by = "file", all = TRUE) }
PER_SAMPLE_COUNTS_merged <- Reduce(f = merge_dt, x = COUNTS)

## Estimate percentage of multiprimer artefacts
PER_SAMPLE_COUNTS_merged[ , 
 PrimerArtefacts_Percent := round(
  PrimerArtefacts_Reads / (PrimerChecked_Reads + PrimerArtefacts_Reads) * 100,
  2)
 ]


### ... update
# .. replace NAs with zero
# .. reorder columns
# .. estimate percentages
# .. add tag-jump summary
# .. add final counts from the Seq table
# .. add positive / negative counts (based on default sample names)

## Prepare per-run stats
PER_RUN_COUNTS_merged <- data.table(
  Total_Number_Of_Reads = RAW$num_seqs,
  Reads_Passed_QC       = QC$num_seqs,
  Reads_Demultiplexed   = sum(PER_SAMPLE_COUNTS_merged$Demultiplexed_Reads, na.rm = TRUE),
  Reads_PrimerChecked   = sum(PER_SAMPLE_COUNTS_merged$PrimerChecked_Reads, na.rm = TRUE)
  )

## Estimate percentages
cat("..Estimating per-run percentages\n")
PER_RUN_COUNTS_merged[ , Percentage_QC := 
  round(Reads_Passed_QC / Total_Number_Of_Reads * 100, 1) ]

PER_RUN_COUNTS_merged[ , Percentage_Demultiplexed := 
  round(Reads_Demultiplexed / Total_Number_Of_Reads * 100, 1) ]

PER_RUN_COUNTS_merged[ , Percentage_Passed := 
  round(Reads_PrimerChecked / Total_Number_Of_Reads * 100, 1) ]


## Export summary stats
write.xlsx(list(
  "per_sample" = PER_SAMPLE_COUNTS_merged,
  "per_run"    = PER_RUN_COUNTS_merged
  ),
  file = "Run_summary.xlsx", colNames = TRUE)


cat("\nAll done.\n")


##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")

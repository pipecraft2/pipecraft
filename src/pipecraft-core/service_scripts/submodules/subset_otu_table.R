#!/usr/bin/env Rscript

## Script to preare a subset of OTU table based on the provided OTU IDs

## Inputs:
# tab-delimited OTU table (columns = OTUs, rows = samples, first column = SampleIDs)
# file with OTU IDs to extract from the OTU table (one ID per row, no header)

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"),  action="store", default=NA, type='character', help="Input OTU table"),
  make_option(c("-s", "--subset"), action="store", default=NA, type='character', help="File with OTU IDs"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Name of the output table")
  )
opt <- parse_args(OptionParser(option_list=option_list))

input_tab <- opt$input
otuids <- opt$subset
output_tab <- opt$output

## Examples
# input_tab <- "DEICODE/full/rclr.tsv"
# otuids <- "OTU_subset_IDs.txt"
# output_tab <- "DEICODE/rclr_subs.tsv"

## Load the table
cat("..Loading OTU table\n")
input_tab <- read.table(file = input_tab,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(input_tab)[1] <- "SampleID"

## Load OTU IDs
cat("..Loading OTU IDs\n")
otuids <- read.table(file = otuids, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(otuids)[1] <- "OTU"

otus <- colnames(input_tab)[ colnames(input_tab) %in% otuids$OTU ]

cat("..Number of OTUs in the OTU table = ", ncol(input_tab) - 1, "\n")
cat("..Number of OTUs in the subset file = ", nrow(otuids), "\n")
cat("..Number of OTUs from the subset file found in the OTU table = ", length(otus), "\n")

## Subset
cat("..Subsetting OTUs\n")
res <- input_tab[, which(colnames(input_tab) %in% c("SampleID", otus)) ]

## Export results
cat("..Exporting OTU table subset\n")
write.table(x = res, file = output_tab,
  quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)

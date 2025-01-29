#!/usr/bin/env Rscript

## Script to perform tag-jump removal
#   by Vladimir Mikryukov

# Input is given as positional arguments:
#   1. OTU table
#   2. f-parameter of UNCROSS 
#   3. p-parameter (e.g., 1.0)

# Outputs:
#  - Tag-jumpfiltered OTU table (`*_TagJumpFilt.txt`)
#  - Plot (`TagJump_plot.pdf`)

# edit 29.01.2025:
#  - add sequences back to the output table (Biostrings and dplyr); needed merge runs process. 

args = commandArgs(trailingOnly = TRUE)

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
theme_set(theme_classic(base_size = 14))

# load env variables
output_dir = Sys.getenv('output_dir')
# Get input file name and create output name
input_file = args[1]
base_name = tools::file_path_sans_ext(basename(input_file))
output_name = paste0(base_name, "_TagJumpFilt.txt")

## Load OTU table
cat(";; Loading input table:", input_file, "\n")
OTUTABW = fread(
  file = input_file,
  sep = "\t", header = TRUE)

## Check table; if 2nd col is sequence, then remove
if (colnames(OTUTABW)[2] == "Sequence") {
  cat(";; 2nd column was 'Sequence', removing for now (adding sequences back later) \n")
  OTUTABW = OTUTABW[, -2]
}

colnames(OTUTABW)[1] = "OTU"

cat(";; Number of ASVs/OTUs: ", nrow(OTUTABW), "\n")
cat(";; Number of samples: ", ncol(OTUTABW) - 1, "\n")

## Convert to long format
cat(";; Converting OTU table to long format\n")
OTUTAB = melt(data = OTUTABW, id.vars = "OTU",
  variable.name = "SampleID", value.name = "Abundance")

## Remove zero-OTUs
OTUTAB = OTUTAB[ Abundance > 0 ]
cat(";; Number of non-zero records: ", nrow(OTUTAB), "\n")

## Estimate total abundance of sequence per plate
cat(";; Estimating total OTU abundance\n")
OTUTAB[ , Total := sum(Abundance, na.rm = TRUE), by = "OTU" ]

## UNCROSS score (with original parameter - take a root from the exp in denominator, to make curves more steep)
uncross_score = function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
  # x = OTU abundance in a sample
  # N = total OTU abundance
  # n = number of samples
  # f = expected cross-talk rate, e.g. 0.01
  # tmin = min score to be considered as cross-talk
  # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

  z = f * N / n               # Expected treshold
  sc = 2 / (1 + exp(x/z)^p)   # t-score
  res = data.table(Score = sc, TagJump = sc >= tmin)
  return(res)
}

## Esimate UNCROSS score
cat(";; Estimating UNCROSS score\n")
OTUTAB = cbind(
  OTUTAB,
  uncross_score(
    x = OTUTAB$Abundance,
    N = OTUTAB$Total,
    n = length(unique(OTUTAB$SampleID)),
    f = as.numeric(args[2]),
    p = as.numeric(args[3])
    )
  )

## Truncate singletons with total OTU abundance > 99 reads
# OTUTAB[ Abundance == 1 & Total > 99  , TagJump := TRUE ]
# OTUTAB[ Abundance == 2 & Total > 999 , TagJump := TRUE ]

cat(";; Number of tag-jumps: ", sum(OTUTAB$TagJump, na.rm = TRUE), "\n")

## Plot
cat(";; Making a plot\n")
PP = ggplot(data = OTUTAB, aes(x = Total, y = Abundance, color = TagJump)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = c("#0C7C59", "#D64933")) +
    labs(x = "Total abundance of OTU, reads", y = "Abundance of OTU in a sample, reads")

cat(";; Exporting a plot\n")
pdf(file = paste0(output_dir, "/TagJump_plot.pdf"), width = 12, height = 9.5, useDingbats = FALSE)
  PP
dev.off()

## TJ stats
cat(";; Calculating tag-jump summary\n")
TJ = data.table(
    Total_reads = sum(OTUTAB$Abundance),
    Number_of_TagJump_Events = sum(OTUTAB$TagJump),
    TagJump_reads = sum(OTUTAB[ TagJump == TRUE ]$Abundance, na.rm = T)
    )

TJ$ReadPercent_removed = with(TJ, (TagJump_reads / Total_reads * 100))

fwrite(x = TJ, file = paste0(output_dir, "/TagJump_stats.txt"), sep = "\t")

## Prepare OTU tables, remove tag-jumps
cat(";; Removing tag-jumps\n")

OTUTAB = OTUTAB[ TagJump == FALSE ]

## Convert to wide format
RES = dcast(data = OTUTAB,
  formula = OTU ~ SampleID,
  value.var = "Abundance", fill = 0)

## Sort rows (by total abundance)
clz = colnames(RES)[-1]
otu_sums = rowSums(RES[, ..clz], na.rm = TRUE)
RES = RES[ order(otu_sums, decreasing = TRUE) ]

## Add sequences back to the table
cat(";; Adding sequences back to the table\n")
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
# read the FASTA file
fasta_sequences = readDNAStringSet(args[4])
sequences = as.character(fasta_sequences)
names(sequences) = names(fasta_sequences)
# add the sequences as 2nd column
RES = RES %>%
  mutate(Sequence = sequences[OTU]) %>%
  select(OTU, Sequence, everything())

## Export table
cat(";; Exporting tag-jump filtered table\n")
cat(paste0(output_dir, output_name))
fwrite(x = RES,
  file = paste0(output_dir, output_name),
  sep = "\t", compress = "none")

cat(";; Done\n")

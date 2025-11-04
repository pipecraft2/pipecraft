#!/usr/bin/env Rscript

## Script to perform tag-jump removal
#   (core function by Vladimir Mikryukov)

# Input is given as positional arguments:
#   1. OTU table (long format; columns - SampleID, OTU, Abundance; with header)
#   2. f-parameter of UNCROSS
#   3. p-parameter (e.g., 1.0)
#   4. fasta file (optional); if provided, sequences will be added back to the output table
#   5. output directory (optional); defaults to $output_dir environment variable

# Outputs:
#  - Tag-jumpfiltered OTU table (`*_TagJumpFilt.txt`)
#  - Plot (`TagJump_plot.pdf`)

args = commandArgs(trailingOnly = TRUE)

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
theme_set(theme_classic(base_size = 14))

# Determine if fasta file is provided (4th arg ends with .fasta or .fa)
has_fasta = FALSE
if (length(args) >= 4) {
    has_fasta = grepl("\\.fa$|\\.fasta$", args[4], ignore.case = TRUE)
    cat(";; has_fasta: ", has_fasta, "\n")
}
# load output directory from args if provided, otherwise from env variable
output_dir = if (length(args) >= (if (has_fasta) 5 else 4)) {
    args[if (has_fasta) 5 else 4]
} else {
    Sys.getenv('output_dir')
}
cat(";; output_dir: ", output_dir, "\n")

# Get input file name and create output name
input_file = args[1]
base_name = tools::file_path_sans_ext(basename(input_file))
output_name = paste0(base_name, "_TagJumpFilt.txt")

## Load OTU table
cat(";; Loading input table:", input_file, "\n")
OTUTAB = fread(
  file = input_file,
  sep = "\t", header = TRUE)

## Check table; if 2nd col is sequence, then remove
if ("Sequence" %in% colnames(OTUTAB)) {
  cat(";; column 'Sequence' is present in the table, removing for now (adding sequences back later) \n")
  OTUTAB[ , Sequence := NULL ]
}

cat(";; Number of ASVs/OTUs: ", uniqueN(OTUTAB$OTU), "\n")
cat(";; Number of samples: ", uniqueN(OTUTAB$SampleID), "\n")

## Remove zero-OTUs
OTUTAB <- OTUTAB[ Abundance > 0 ]
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
    n = uniqueN(OTUTAB$SampleID),
    f = as.numeric(args[2]),
    p = as.numeric(args[3])
    )
  )

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

## Remove working columns
OTUTAB[ , c("Total", "Score", "TagJump") := NULL ]

## Sort rows (by sample and abundance)
setorder(OTUTAB, SampleID, -Abundance, OTU)

## Add sequences back to the table (if FASTA file provided)
if (has_fasta) {
    cat(";; Adding sequences back to the table\n")
    suppressMessages(library(Biostrings))
    suppressMessages(library(dplyr))
    # read the FASTA file
    cat(";; Loading input FASTA:", args[4], "\n")
    fasta_sequences = readDNAStringSet(args[4])
    sequences = as.character(fasta_sequences)
    names(sequences) = names(fasta_sequences)
    # add the sequences as 2nd column
    RES = RES %>%
    mutate(Sequence = sequences[OTU]) %>%
    select(OTU, Sequence, everything())
}

## Export table
cat(";; Exporting tag-jump filtered table", 
    if(length(args) >= 4) "(2nd col is sequence)" else "", "\n")
cat(paste0(output_dir, "/", output_name))
fwrite(x = RES, file = paste0(output_dir, "/", output_name), sep = "\t", compress = "none")

cat(";; Done\n")
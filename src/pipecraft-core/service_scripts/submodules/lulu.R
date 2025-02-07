#!/usr/bin/env Rscript
library(devtools)

# load variables
min_ratio_type = Sys.getenv('min_ratio_type')
min_ratio = as.numeric(Sys.getenv('min_ratio'))
min_match = as.numeric(Sys.getenv('min_match'))
min_rel_cooccurence = as.numeric(Sys.getenv('min_rel_cooccurence'))
fasta_file = Sys.getenv('input_fasta')
otu_table = Sys.getenv('otu_table')

# load OTU table
cat(";; Loading OTU table\n")
otutable = read.table(otu_table, 
					header = T, 
					row.names = 1, 
					sep = "\t", 
					check.names = FALSE)
# load match list
cat(";; Loading match list\n")
matchlist_name = read.table(file.path("/input/lulu_out", 
	"match_list.lulu"))

# check for "Sequence" column in OTU table
	# remove for LULU, but insert it back in the end based on fasta file
if (colnames(otutable)[1] == "Sequence") {
  cat(";; 2nd column was 'Sequence', removing for now (adding sequences back later) \n")
  otutable = otutable[, -1]
}

# run lulu
cat(";; ")
cat(";; Running LULU with the following parameters:\n")
cat(";; minimum_ratio_type =", min_ratio_type, "\n")
cat(";; minimum_ratio =", min_ratio, "\n")
cat(";; minimum_match =", min_match, "\n")
cat(";; minimum_relative_cooccurence =", min_rel_cooccurence, "\n")
cat(";; input otu table =", otu_table, "\n")
cat(";; ")
curated_result <- lulu::lulu(otutable, matchlist_name, 
	minimum_ratio_type = min_ratio_type, 
	minimum_ratio = min_ratio, 
	minimum_match = min_match, 
	minimum_relative_cooccurence = min_rel_cooccurence)
cat(";; ")
cat(";; LULU completed;;\n")
# output OTU table
output_otutable = curated_result$curated_table

# get OTU ids from curated table 
cat(";; Outputting OTU ids\n")
write.table(rownames(output_otutable), file ="/input/lulu_out/lulu_out_OTUids.txt", row.names = FALSE, quote = FALSE)

## Add sequences back to the table
cat(";; Adding sequences back to the table\n")
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
# read the FASTA file
cat(";; Loading input FASTA:", fasta_file, "\n")
fasta_sequences = readDNAStringSet(fasta_file)
sequences = as.character(fasta_sequences)
names(sequences) = names(fasta_sequences)
# add the sequences as 2nd column
output_otutable = output_otutable %>%
	mutate(Sequence = sequences[row.names(output_otutable)]) %>%
	select(Sequence, everything())

# get the basename of the input otu table
otu_table_basename = sub("\\.[^.]*$", "", basename(otu_table))

#write post-clustered OTU table to file
cat(";; Writing post-clustered OTU table to file\n")
write.table(output_otutable, file = file.path("/input/lulu_out", paste0(otu_table_basename, ".lulu.txt")), 
	sep = "\t", 
	col.names = NA, 
	row.names = TRUE, 
	quote = FALSE)

cat(";; Writing discarded units to file\n")
write.table(curated_result$discarded_otus, file ="/input/lulu_out/discarded_units.lulu", col.names = FALSE, quote = FALSE)

# Remove original OTU table from $outputdir
if (file.exists((file.path("/input/lulu_out", "OTU_tab_for_lulu.txt")))) {
    unlink((file.path("/input/lulu_out", "OTU_tab_for_lulu.txt")))
}

cat(";; Done \n")
#DONE, proceed with lulu.sh to clean up make readme

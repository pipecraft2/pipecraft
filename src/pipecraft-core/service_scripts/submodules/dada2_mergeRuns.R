#!/usr/bin/env Rscript

# DADA2 merge sequencing runs

# load dada2
library('dada2')
# print DADA2 version
cat(";; DADA2, version", base::toString(packageVersion("dada2")), "\n")

# get input tables (txt files)
output_feature_table=Sys.getenv("output_feature_table")
cat(";; output_feature_table: ", output_feature_table, "\n")
output_fasta=Sys.getenv("output_fasta")
output_dir=Sys.getenv("output_dir")

# Check if output_feature_table is empty
if (output_feature_table == "") {
    stop("Error: No input tables provided in output_feature_table environment variable")
}

# Split the comma-separated string into a vector of file paths
input_tables = strsplit(output_feature_table, ",")[[1]]
input_tables = trimws(input_tables)  # Remove any whitespace

# Check if files exist
missing_files = input_tables[!file.exists(input_tables)]
if (length(missing_files) > 0) {
    stop("Error: Following files not found: ", paste(missing_files, collapse=", "))
}

# Initialize empty list for tables
input_tables_list = list()
# Read each table and store in the list
for (i in 1:length(input_tables)) { # nolint
    input_tables_list[[i]] = read.table(input_tables[i], 
                                  header = T, 
                                  row.names = 1, 
                                  check.names = F, 
                                  sep = "\t")
}

# Format the txt ASV table for dada2 mergeSequenceTables
for (i in seq_along(input_tables)) {
    print(paste0(";; Formatting table ", i, "."))
    tab = input_tables_list[[i]]
    # sequences are row names
    rownames(tab) = tab$Sequence
    # remove col "Sequence"
    tab = tab[, !(colnames(tab) == "Sequence")]
    #transpose the tab
    tab = t(tab)
    # store modified table back in list
    input_tables_list[[i]] = tab
}


# MERGE RUNS with mergeSequenceTables
cat(";; Merging tables...\n")
merged_table = mergeSequenceTables(tables=input_tables_list, repeats = "sum")

###format and save ASV table and ASVs.fasta
#sequence headers
cat(";; Formatting ASV table and ASVs.fasta...\n")
asv_seqs = colnames(merged_table)
asv_headers = openssl::sha1(asv_seqs)

#transpose sequence table
tmerged_table = t(merged_table)
#add sequences to 1st column
tmerged_table = cbind(row.names(tmerged_table), tmerged_table)
colnames(tmerged_table)[1] = "Sequence"
#row names as sequence headers
row.names(tmerged_table) = asv_headers

#write ASVs.fasta to path_ASVs
cat(";; Writing ASVs.fasta...\n")
asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
write(asv_fasta, file.path(output_dir, "ASVs.fasta"))
#write ASVs table to path_ASVs
cat(";; Writing ASVs table...\n")
write.table(tmerged_table, file.path(output_dir, "ASVs_table.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

cat(";; Done!\n")


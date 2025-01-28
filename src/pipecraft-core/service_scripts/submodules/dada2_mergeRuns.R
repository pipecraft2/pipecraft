#!/usr/bin/env Rscript

# DADA2 merge sequencing runs

# load dada2
library('dada2')

# print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

# get input tables (txt files)
input_tables = list()
tab.1 = read.table(input_table1, 
            header = T, row.names = 1, check.names = F, sep = "\t")


# Format the txt ASV table for mergeSequenceTables
i=1
for (i in 1:length(input_tables)){
    print(i)
    tab = get(paste0("tab.", i))
    # sequences are row names
    rownames(tab) = tab$Sequence
    # remove col "Sequence"
    tab = tab[, !(colnames(tab) == "Sequence")]
    #transpose the tab
    tab = t(tab)
    print(rownames(tab))
    # assign the modified table back to its original name
    assign(paste0("tab.", i), tab)
}



# MERGE RUNS
merged_table = mergeSequenceTables(tables=input_tables, repeats = "sum")

###format and save ASV table and ASVs.fasta
#sequence headers
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
asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
write(asv_fasta, file.path(path_ASVs, "ASVs.fasta"))
#write ASVs table to path_ASVs
write.table(tmerged_table, file.path(path_ASVs, "ASVs_table.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)



# Remove chimeras (cuz those are not removed from the PacBio data)
print("Removing chimeras...")
ASV_tab.nochim = removeBimeraDenovo(merged_table, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(ASV_tab.nochim, "ASV_tab.nochim.rds")

###format and save ASV table and ASVs.fasta
#sequence headers
asv_seqs = colnames(ASV_tab.nochim)
asv_headers = openssl::sha1(asv_seqs)

#transpose sequence table
tASV_tab.nochim = t(ASV_tab.nochim)
#add sequences to 1st column
tASV_tab.nochim = cbind(row.names(tASV_tab.nochim), tASV_tab.nochim)
colnames(tASV_tab.nochim)[1] = "Sequence"
#row names as sequence headers
row.names(tASV_tab.nochim) = asv_headers

#write ASVs.fasta to path_ASVs
asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
write(asv_fasta, file.path(path_ASVs, "ASVs.nochim.fasta"))
#write ASVs table to path_ASVs
write.table(tASV_tab.nochim, file.path(path_ASVs, "ASVs_table.nochim.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)


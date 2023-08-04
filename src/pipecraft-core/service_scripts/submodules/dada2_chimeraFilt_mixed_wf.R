#!/usr/bin/env Rscript

#DADA2 chimera filtering for DADA2 full workflow. In this step, also ASV table is generated (and ASVs.fasta)

#load dada2
library("dada2")
library("base")

#load env variables
readType = Sys.getenv('readType')
fileFormat = Sys.getenv('fileFormat')
dataFormat = Sys.getenv('dataFormat')
workingDir = Sys.getenv('workingDir')

#load variables
method = Sys.getenv('method')

#check for output dir and delete if needed
if (dir.exists("/input/chimeraFiltered_out.dada2")) {
    unlink("/input/chimeraFiltered_out.dada2", recursive=TRUE)
}
if (dir.exists("/input/ASVs_out.dada2")) {
    unlink("/input/ASVs_out.dada2", recursive=TRUE)
}
#create output dir
path_results = "/input/chimeraFiltered_out.dada2"
dir.create(path_results)
path_ASVs = "/input/ASVs_out.dada2"
dir.create(path_ASVs)


### combine fwd_orient and rev_orient runs with mergeSequenceTables
cat("| Combining fwd_orient and rev_orient runs with mergeSequenceTables | ")
#load fwd and rev tables
fwd = readRDS("/input/primersCut_out/fwd_orient/denoised_assembled.dada2/ASVs_table.denoised.rds")
rev = readRDS("/input/primersCut_out/rev_orient/denoised_assembled.dada2/ASVs_table.denoised.rds")
### RevComp ASV sequences in rev tables 
colnames(rev) = dada2:::rc(colnames(rev))

# Merge tables
    # repeats = "sum" -> samples with the same name are summed together in the merged table.
ASV_tab = mergeSequenceTables(fwd, rev, repeats = "sum") 
saveRDS(ASV_tab, file.path(path_ASVs, "ASVs_table.denoised.rds"))    

#remove chimeras
cat("| Removing chimeras with removeBimeraDenovo | ")
ASV_tab.nochim <- removeBimeraDenovo(ASV_tab, method = method, multithread = FALSE, verbose = TRUE)

#save rds
saveRDS(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.denoised.nochim.rds"))

#seq count summary
seq_count <- cbind(rowSums(ASV_tab), rowSums(ASV_tab.nochim))
colnames(seq_count) <- c("input(merged)", "chimeraFiltered")
write.table(seq_count, file.path(path_results, "seq_count_summary.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

###format and save ASVs_table.txt and ASVs.fasta
#sequence headers
asv_seqs = colnames(ASV_tab.nochim)
asv_size = colSums(ASV_tab.nochim)
asv_headers = openssl::sha1(asv_seqs)

#transpose sequence table
ASV_tab.nochim = t(ASV_tab.nochim)
#add sequences to 1st column
ASV_tab.nochim = cbind(row.names(ASV_tab.nochim), ASV_tab.nochim)
colnames(ASV_tab.nochim)[1] = "Sequence"
#row names as sequence headers
row.names(ASV_tab.nochim) = asv_headers

#write ASVs.fasta to path_ASVs
asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
write(asv_fasta, file.path(path_ASVs, "ASVs.fasta"))
#write ASVs table to path_ASVs
write.table(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#Loop through each sample in the table and write per-sample fasta files containin non-chimeric ASVs (start with 2nd col [1st is Sequence])
for (i in 2:length(colnames(ASV_tab.nochim))){
    sample_name = colnames(ASV_tab.nochim)[i]
    sample_file = paste(sample_name, "chimFilt_ASVs.fasta", sep = ".") 
    j = 0
    for (abundance in ASV_tab.nochim[,i]){
        j = j + 1
        if (abundance != 0){
            #seq header and abundance
            header = paste(">", row.names(ASV_tab.nochim)[j], ";size=", abundance, sep = "")
            write(header, file.path(path_results, sample_file), append = TRUE)
            #sequence
            seq = ASV_tab.nochim[j, 1]
            write(seq, file.path(path_results, sample_file), append = TRUE)
        }
    }
}

#DONE, proceed with chimera_filtering_dada2_wf.sh to clean up make readme

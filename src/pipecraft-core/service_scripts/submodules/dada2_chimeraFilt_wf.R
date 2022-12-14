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

#load data
ASV_tab = readRDS(file.path(workingDir, "ASVs_table.denoised.rds"))

#remove chimeras
ASV_tab.nochim <- removeBimeraDenovo(ASV_tab, method = method, multithread = FALSE, verbose = FALSE)

#save rds
saveRDS(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.denoised.nochim.rds"))

#seq count summary
qfilt = readRDS("/input/qualFiltered_out/quality_filtered.rds")
sample_names = readRDS("/input/qualFiltered_out/sample_names.rds")

no_of_ASVs_list = list() #add ASVs per sample count
for (i in 1:nrow(ASV_tab.nochim)){
    no_of_ASVs = sum(ASV_tab.nochim[i,] > 0)
    no_of_ASVs_list = append(no_of_ASVs_list, no_of_ASVs, after = length(no_of_ASVs_list))
}
ASV_counts = data.frame(no_of_ASVs_list, check.names = FALSE, row.names = "")
colnames(ASV_counts) = sample_names
seq_count <- cbind(rowSums(ASV_tab), rowSums(ASV_tab.nochim), t(ASV_counts))
colnames(seq_count) <- c("input(merged)", "chimeraFiltered", "no.of ASVs")
rownames(seq_count) <- sample_names
write.table(seq_count, file.path(path_results, "seq_count_summary.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

###format and save ASVs_table.txt and ASVs.fasta
#sequence headers
asv_seqs = colnames(ASV_tab.nochim)
asv_size = colSums(ASV_tab.nochim)
asv_headers = vector(dim(ASV_tab.nochim)[2], mode="character")
for (i in 1:dim(ASV_tab.nochim)[2]) {
    asv_headers[i] = paste(">ASV", i, ";size=", asv_size[i], sep="")
}
#transpose sequence table
ASV_tab.nochim = t(ASV_tab.nochim)
#add sequences to 1st column
ASV_tab.nochim = cbind(row.names(ASV_tab.nochim), ASV_tab.nochim)
colnames(ASV_tab.nochim)[1] = "Sequence"
#row names as sequence headers
row.names(ASV_tab.nochim) = sub(">", "", asv_headers)
#write ASVs.fasta to path_ASVs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(path_ASVs, "ASVs.fasta"))
#write ASVs table to path_ASVs
write.table(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#Loop through each sample in the table and write per-sample fasta files containin non-chimeric ASVs
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

#remove not needed .rds objects
if (file.exists("/input/qualFiltered_out/filtFs.rds")) {
    unlink("/input/qualFiltered_out/*.rds")
}

#DONE, proceed with chimera_filtering_dada2_wf.sh to clean up make readme

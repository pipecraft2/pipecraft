#!/usr/bin/env Rscript

#DADA2 chimera filtering for DADA2 full workflow. In this step, also ASV table is generated (and ASVs.fasta)

#load dada2
library("dada2")
library("base")

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

#load env variables
readType = Sys.getenv("readType")
fileFormat = Sys.getenv("fileFormat")
dataFormat = Sys.getenv("dataFormat")
workingDir = Sys.getenv("workingDir")
cat(";; workingDir = ", workingDir, "\n")

#load variables
method = Sys.getenv("method")
allowOneOff = as.logical(Sys.getenv("allowOneOff"))
minOneOffParentDistance = as.numeric(Sys.getenv("minOneOffParentDistance"))
maxShift = as.numeric(Sys.getenv("maxShift"))
path_results = Sys.getenv("output_dir1")
path_ASVs = Sys.getenv("output_dir2")

cat(";; Settings:\n")
cat(";; method = ", method, "\n")
cat(";; allowOneOff = ", allowOneOff, "\n")
cat(";; minOneOffParentDistance = ", minOneOffParentDistance, "\n")
cat(";; maxShift = ", maxShift, ";; \n")

#check for output dirs and delete if needed
if (dir.exists(path_results)) {
    unlink(path_results, recursive = TRUE)
}
dir.create(path_results)
if (dir.exists(path_ASVs)) {
    unlink(path_ASVs, recursive = TRUE)
}
dir.create(path_ASVs)

#load data
ASV_tab = readRDS(file.path(workingDir, "ASVs_table.denoised.rds"))

#remove chimeras
cat(";; Removing chimeras with removeBimeraDenovo function ;; ")
ASV_tab.nochim = removeBimeraDenovo(ASV_tab, method = method, multithread = TRUE, verbose = TRUE)

#save rds
saveRDS(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.denoised.nochim.rds"))
saveRDS(ASV_tab, file.path(path_ASVs, "ASVs_table.denoised.rds"))

#seq count summary
    # if sample_names.rds is not in the workingDir, copy it from the qualFiltered_out dir
if (!file.exists(file.path(workingDir, "sample_names.rds"))) {
    qualFiltered_out_dir = file.path(dirname(workingDir), "qualFiltered_out")
    file.copy(from = file.path(qualFiltered_out_dir, "sample_names.rds"),
              to = file.path(workingDir, "sample_names.rds"))
}
sample_names = readRDS(file.path(workingDir, "sample_names.rds"))
cat(";; sample names = ", sample_names)

no_of_ASVs_list = list() #add ASVs per sample count
for (i in 1:nrow(ASV_tab.nochim)){
    no_of_ASVs = sum(ASV_tab.nochim[i,] > 0)
    no_of_ASVs_list = append(no_of_ASVs_list, no_of_ASVs, after = length(no_of_ASVs_list))
}
ASV_counts = data.frame(no_of_ASVs_list, check.names = FALSE, row.names = "")
colnames(ASV_counts) = sample_names
seq_count <- cbind(rowSums(ASV_tab), rowSums(ASV_tab.nochim), t(ASV_counts))
colnames(seq_count) <- c("input(denoised)", "chimeraFiltered", "no.of ASVs")
rownames(seq_count) <- sample_names
write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), row.names = TRUE, quote = FALSE)

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
asv_fasta <- c(rbind(paste(">", asv_headers, sep = ""), asv_seqs))
write(asv_fasta, file.path(path_ASVs, "ASVs.fasta"))
#write ASVs table to path_ASVs
write.table(ASV_tab.nochim, file.path(path_ASVs, "ASVs_table.txt"),
            sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#Loop through each sample in the table and write per-sample fasta files containin non-chimeric ASVs
cat(";; Writing per-sample fasta files containin non-chimeric ASVs ")
for (i in 2:length(colnames(ASV_tab.nochim))){
    sample_name = colnames(ASV_tab.nochim)[i]
    sample_file = paste(sample_name, "chimFilt_ASVs.fasta", sep = ".")
    j = 0
    for (abundance in ASV_tab.nochim[, i]){
        j = j + 1
        if (abundance != 0) {
            #seq header and abundance
            header = paste(">", row.names(ASV_tab.nochim)[j], ";size=", abundance, sep = "")
            write(header, file.path(path_results, sample_file), append = TRUE)
            #sequence
            seq = ASV_tab.nochim[j, 1]
            write(seq, file.path(path_results, sample_file), append = TRUE)
        }
    }
}
cat(";; DONE ")
#DONE, proceed with chimera_filtering_dada2_wf.sh to clean up make readme

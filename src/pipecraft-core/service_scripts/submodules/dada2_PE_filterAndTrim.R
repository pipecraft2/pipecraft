#!/usr/bin/env Rscript

#DADA2 PE data quality filtering.

#load dada2
library('dada2')

#load env variables
fileFormat = Sys.getenv('fileFormat')

#load variables
read_R1 = Sys.getenv('read_R1')
read_R2 = Sys.getenv('read_R2')
samp_ID = Sys.getenv('samp_ID')
maxEE = as.numeric(Sys.getenv('maxEE'))
maxN = as.numeric(Sys.getenv('maxN'))
truncQ = as.numeric(Sys.getenv('truncQ'))
truncLen_R1 = as.numeric(Sys.getenv('truncLen'))
truncLen_R2 = as.numeric(Sys.getenv('truncLen_R2'))
minLen = as.numeric(Sys.getenv('minLen'))
maxLen = as.numeric(Sys.getenv('maxLen'))
minQ = as.numeric(Sys.getenv('minQ'))

#output path
path_results = "/input/qualFiltered_out"

#define input and output file paths
fnFs = sort(list.files(pattern = read_R1, full.names = TRUE))
fnRs = sort(list.files(pattern = read_R2, full.names = TRUE))
print(fnFs)
print(fnRs)
#sample names
sample_names = sapply(strsplit(basename(fnFs), samp_ID), `[`, 1)

#filtered files path
filtFs = file.path(path_results, paste0(sample_names, "_R1_filt.", "fastq"))
filtRs = file.path(path_results, paste0(sample_names, "_R2_filt.", "fastq"))
names(filtFs) = sample_names
names(filtRs) = sample_names
print(filtFs)
print(filtRs)

#quality filter
qfilt = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    maxN = maxN, 
                    maxEE = c(maxEE, maxEE), 
                    truncQ = truncQ,  
                    truncLen = c(truncLen_R1, truncLen_R2),
                    maxLen = maxLen, 
                    minLen = minLen, 
                    minQ = minQ, 
                    rm.phix = TRUE, 
                    matchIDs = FALSE,
                    compress = FALSE, 
                    multithread = TRUE)
saveRDS(qfilt, file.path(path_results, "quality_filtered.rds"))

#seq count summary
getN <- function(x) sum(getUniques(x))
seq_count <- cbind(qfilt)
colnames(seq_count) <- c("input", "qualFiltered")
rownames(seq_count) <- sample_names
write.table(seq_count, file.path(path_results, "seq_count_summary.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#save R objects for assembly process
R1qf = sort(list.files(path_results, pattern = "_R1_filt.", full.names = TRUE))
R2qf = sort(list.files(path_results, pattern = "_R2_filt.", full.names = TRUE))
sample_names = sapply(strsplit(basename(R1qf), "_R1_filt."), `[`, 1)
saveRDS(R1qf, file.path(path_results, "filtFs.rds"))
saveRDS(R2qf, file.path(path_results, "filtRs.rds"))
saveRDS(sample_names, file.path(path_results, "sample_names.rds"))

#DONE, proceed with quality_filtering_paired_end_dada2.sh to clean up make readme

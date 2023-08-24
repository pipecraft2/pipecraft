#!/usr/bin/env Rscript

#DADA2 PE data quality filtering.

#load dada2
library('dada2')

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

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

#check if gz files are provided; if yes then produce also gz compressed files.
is_gz = strsplit(fileFormat, split="\\.")[[1]][-1]
if (identical(is_gz, character(0)) != "TRUE") {
    if (is_gz == "gz") {
        compress = TRUE
    } else {
        compress = FALSE
    }
} else {
    compress = FALSE
}

#output path
path_results = Sys.getenv('output_dir')

#define input and output file paths
fnFs = sort(list.files(pattern = read_R1, full.names = TRUE))
fnRs = sort(list.files(pattern = read_R2, full.names = TRUE))
#sample names
sample_names = sapply(strsplit(basename(fnFs), samp_ID), `[`, 1)
cat(";; sample names = ", sample_names, "\n")

#filtered files path
filtFs = file.path(path_results, paste0(sample_names, "_R1.", fileFormat))
filtRs = file.path(path_results, paste0(sample_names, "_R2.", fileFormat))
names(filtFs) = sample_names
names(filtRs) = sample_names

#quality filter
cat(";; Quality filtering with filterAndTrim function")
qfilt = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    maxN = maxN, 
                    maxEE = c(maxEE, maxEE), 
                    truncQ = truncQ,  
                    truncLen = c(truncLen_R1, truncLen_R2),
                    maxLen = maxLen, 
                    minLen = minLen, 
                    minQ = minQ, 
                    rm.phix = TRUE, 
                    matchIDs = "FALSE",
                    compress = compress, 
                    multithread = TRUE)
saveRDS(qfilt, file.path(path_results, "quality_filtered.rds"))

#seq count summary
getN <- function(x) sum(getUniques(x))
seq_count <- cbind(qfilt)
colnames(seq_count) <- c("input", "qualFiltered")
rownames(seq_count) <- sample_names
write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), row.names = TRUE, quote = FALSE)

#save R objects for assembly process
R1qf = sort(list.files(path_results, pattern = "_R1.", full.names = TRUE))
R2qf = sort(list.files(path_results, pattern = "_R2.", full.names = TRUE))
sample_names = sapply(strsplit(basename(R1qf), "_R1."), `[`, 1)
saveRDS(R1qf, file.path(path_results, "filtFs.rds"))
saveRDS(R2qf, file.path(path_results, "filtRs.rds"))
saveRDS(sample_names, file.path(path_results, "sample_names.rds"))
cat(";; DONE ")
#DONE, proceed with quality_filtering_paired_end_dada2.sh to clean up make readme

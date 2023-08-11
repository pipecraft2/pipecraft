#!/usr/bin/env Rscript

#DADA2 SE data quality filtering.

#load dada2
library('dada2')

#load env variables
fileFormat = Sys.getenv('fileFormat')
#output path
path_results = Sys.getenv('output_dir')
#load variables
maxEE = as.numeric(Sys.getenv('maxEE'))
maxN = as.numeric(Sys.getenv('maxN'))
truncQ = as.numeric(Sys.getenv('truncQ'))
truncLen_R1 = as.numeric(Sys.getenv('truncLen'))
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

#define input and output file paths
fnFs = sort(list.files(pattern = fileFormat, full.names = TRUE))
#sample names
sample_names = sapply(strsplit(basename(fnFs), paste0(".", fileFormat)), `[`, 1)
cat("sample names = ", sample_names, "\n")

#filtered files path
qFiltered = file.path(path_results, paste0(sample_names, ".", fileFormat))
names(qFiltered) = sample_names

#quality filter
qfilt = filterAndTrim(fnFs, qFiltered, 
                    maxN = maxN, 
                    maxEE = maxEE, 
                    truncQ = truncQ,  
                    truncLen = truncLen_R1,
                    maxLen = maxLen, 
                    minLen = minLen, 
                    minQ = minQ, 
                    rm.phix = TRUE, 
                    compress = compress, 
                    multithread = TRUE, verbose = TRUE)
saveRDS(qfilt, file.path(path_results, "quality_filtered_read_count.rds"))

#seq count summary
getN <- function(x) sum(getUniques(x))
seq_count <- cbind(qfilt)
colnames(seq_count) <- c("input", "qualFiltered")
rownames(seq_count) <- sample_names
write.csv(seq_count, file.path(path_results, "seq_count_summary.txt"), row.names = TRUE, quote = FALSE)

#save R objects for denoising
filtered = sort(list.files(path_results, pattern = fileFormat, full.names = TRUE))
sample_names = sapply(strsplit(basename(filtered), paste0(".", fileFormat)), `[`, 1)
saveRDS(filtered, file.path(path_results, "qFiltered.rds"))
saveRDS(sample_names, file.path(path_results, "sample_names.rds"))

#DONE, proceed with quality_filtering_single_end_dada2.sh to clean up make readme

#!/usr/bin/env Rscript

#DADA2 denoising and merging paired-end data

#load dada2
library('dada2')

#load env variables
fileFormat = Sys.getenv('fileFormat')

#load  variables
read_R1 = Sys.getenv('read_R1')
read_R2 = Sys.getenv('read_R2')
samp_ID = Sys.getenv('samp_ID')
minOverlap = as.numeric(Sys.getenv('minOverlap'))
maxMismatch = as.numeric(Sys.getenv('maxMismatch'))
trimOverhang = Sys.getenv('trimOverhang')
justConcatenate = Sys.getenv('justConcatenate')
pool = Sys.getenv('pool')
selfConsist = Sys.getenv('selfConsist')
qualityType = Sys.getenv('qualityType')

#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (pool == "false" || pool == "FALSE"){
    pool = FALSE
}
if (pool == "true" || pool == "TRUE"){
    pool = TRUE
}
if (selfConsist == "false" || selfConsist == "FALSE"){
    selfConsist = FALSE
}
if (selfConsist == "true" || selfConsist == "TRUE"){
    selfConsist = TRUE
}
if (trimOverhang == "false" || trimOverhang == "FALSE"){
    trimOverhang = FALSE
}
if (justConcatenate == "false" || justConcatenate == "FALSE"){
    justConcatenate = FALSE
}
if (trimOverhang == "true" || trimOverhang == "TRUE"){
    trimOverhang = TRUE
}
if (justConcatenate == "true" || justConcatenate == "TRUE"){
    justConcatenate = TRUE
}

#output path
path_results = "/input/denoised_assembled.dada2"

#define input file paths
fnFs = sort(list.files(pattern = read_R1, full.names = TRUE))
fnRs = sort(list.files(pattern = read_R2, full.names = TRUE))
#sample names
sample_names = sapply(strsplit(basename(fnFs), samp_ID), `[`, 1)
cat("sample names = ", sample_names, "\n")

cat("#Performing DADA2 denoising;")

#Learn the error rates
errF = learnErrors(fnFs, multithread = FALSE)
errR = learnErrors(fnRs, multithread = FALSE)

#Error rate figures
pdf(file.path(path_results, "Error_rates_R1.pdf"))
    print( plotErrors(errF) )
dev.off()
pdf(file.path(path_results, "Error_rates_R2.pdf"))
    print( plotErrors(errR) )
dev.off()

#dereplicate
derepFs = derepFastq(fnFs, qualityType = qualityType)
derepRs = derepFastq(fnRs, qualityType = qualityType)
print("derepFastq DONE")

#denoise
dadaFs = dada(derepFs, err = errF, pool = pool, selfConsist = selfConsist, multithread = FALSE)
dadaRs = dada(derepRs, err = errR, pool = pool, selfConsist = selfConsist, multithread = FALSE)
print("denoising DONE")

#merge paired-end reads
cat("#Merging data with mergePairs;")
merge = mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                            maxMismatch = maxMismatch,
                            minOverlap = minOverlap,
                            justConcatenate = justConcatenate,
                            trimOverhang = trimOverhang)

### WRITE temporary PER-SAMPLE DENOISED and MERGED ASV FASTA FILES
#make sequence table
ASV_tab = makeSequenceTable(merge)
#write RDS object
#saveRDS(ASV_tab, (file.path(path_results, "ASVs_table.denoised.rds")))

#sequence headers
asv_seqs = colnames(ASV_tab)
asv_headers = openssl::sha1(asv_seqs) #header as sha1

#transpose sequence table
ASV_tab = t(ASV_tab)
#add sequences to 1st column
ASV_tab = cbind(row.names(ASV_tab), ASV_tab)
colnames(ASV_tab)[1] = "Sequence"
#row names as sequence headers
row.names(ASV_tab) = asv_headers

#write ASVs table
#write.table(ASV_tab, file.path(path_results, "ASVs_table.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#Loop through each sample in the table and write per-sample fasta files
for (i in 2:length(colnames(ASV_tab))){
    sample_name = colnames(ASV_tab)[i]
    sample_file = paste(sample_name, "ASVs.fasta", sep = ".") 
    j = 0
    for (abundance in ASV_tab[,i]){
        j = j + 1
        if (abundance != 0){
            #seq header and abundance
            header = paste(">", row.names(ASV_tab)[j], ";size=", abundance, sep = "")
            write(header, file.path(path_results, sample_file), append = TRUE)
            #sequence
            seq = ASV_tab[j, 1]
            write(seq, file.path(path_results, sample_file), append = TRUE)
        }
    }
}

#seq count summary
getN <- function(x) sum(getUniques(x))
seq_count <- cbind(sapply(derepFs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merge, getN))
colnames(seq_count) <- c("input", "denoised_R1", "denoised_R2", "merged")
rownames(seq_count) <- sample_names
write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), row.names = TRUE, quote = FALSE)
print("DONE")
#DONE, proceed with assemble_paired_end_dada2.sh to clean up make readme

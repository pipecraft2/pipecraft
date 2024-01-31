#!/usr/bin/env Rscript

#DADA2 denoising and merging paired-end data

#load dada2
library('dada2')

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

#load env variables
fileFormat = Sys.getenv('fileFormat')

#load  variables
read_R1 = Sys.getenv('read_R1')
read_R2 = Sys.getenv('read_R2')
minOverlap = as.numeric(Sys.getenv('minOverlap'))
maxMismatch = as.numeric(Sys.getenv('maxMismatch'))
trimOverhang = Sys.getenv('trimOverhang')
justConcatenate = Sys.getenv('justConcatenate')
pool = Sys.getenv('pool')
qualityType = Sys.getenv('qualityType')
band_size = as.numeric(Sys.getenv('BAND_SIZE'))
#setDadaOpt() settings
omegaa = as.numeric(Sys.getenv('OMEGA_A'))
omegap = as.numeric(Sys.getenv('OMEGA_P'))
omegac= as.numeric(Sys.getenv('OMEGA_C'))
detect_singletons = Sys.getenv('DETECT_SINGLETONS')
band_size = as.numeric(Sys.getenv('BAND_SIZE'))


#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (pool == "false" || pool == "FALSE"){
    pool = FALSE
}
if (pool == "true" || pool == "TRUE"){
    pool = TRUE
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
if (detect_singletons == "false" || detect_singletons == "FALSE"){
    detect_singletons = FALSE
}
if (detect_singletons == "true" || detect_singletons == "TRUE"){
    detect_singletons = TRUE
}

#Set DADA options
setDadaOpt(OMEGA_A = omegaa, OMEGA_P = omegap, OMEGA_C = omegac, DETECT_SINGLETONS = detect_singletons, BAND_SIZE = band_size)
cat(";; BAND_SIZE = ", band_size, "\n")  

#output path
path_results = "/input/denoised_assembled.dada2"

#define input file paths
fnFs = sort(list.files(pattern = read_R1, full.names = TRUE))
fnRs = sort(list.files(pattern = read_R2, full.names = TRUE))
#sample names
sample_names = sapply(strsplit(basename(fnFs), read_R1), `[`, 1)
cat(";; sample names = ", sample_names, "\n")

#Learn the error rates
cat(";; Performing DADA2 denoising ;; ")
errF = learnErrors(fnFs, multithread = TRUE)
errR = learnErrors(fnRs, multithread = TRUE)

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
dadaFs = dada(derepFs, err = errF, pool = pool, multithread = TRUE)
dadaRs = dada(derepRs, err = errR, pool = pool, multithread = TRUE)
print(";; denoising DONE")

#merge paired-end reads
cat(";; Merging data with mergePairs function ")
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
cat(";; Writing per-sample fasta files ")
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
cat(";; DONE")
#DONE, proceed with assemble_paired_end_dada2.sh to clean up make readme

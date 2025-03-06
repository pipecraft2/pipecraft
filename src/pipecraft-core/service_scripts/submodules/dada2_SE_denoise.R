#!/usr/bin/env Rscript

#DADA2 denoising single-end data.

#load dada2
library('dada2')

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

#load env variables
readType = Sys.getenv('readType')
fileFormat= Sys.getenv('fileFormat')
dataFormat = Sys.getenv('dataFormat')
workingDir = Sys.getenv('workingDir')

#load  variables
errorEstFun = Sys.getenv('errorEstFun')
pool = Sys.getenv('pool')
qualityType = Sys.getenv('qualityType')
#setDadaOpt() settings
omegaa = as.numeric(Sys.getenv('OMEGA_A'))
omegap = as.numeric(Sys.getenv('OMEGA_P'))
omegac= as.numeric(Sys.getenv('OMEGA_C'))
detect_singletons = Sys.getenv('DETECT_SINGLETONS')
band_size = as.numeric(Sys.getenv('BAND_SIZE'))


cat(";; errorEstimationFunction = ", errorEstFun, "\n")
cat(";; BAND_SIZE = ", band_size, "\n")

#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (pool == "false" || pool == "FALSE"){
    pool = FALSE
}
if (pool == "true" || pool == "TRUE"){
    pool = TRUE
}
if (errorEstFun == "PacBioErrfun"){
    errorEstFun = PacBioErrfun
} else {
    errorEstFun = loessErrfun
}
if (detect_singletons == "false" || detect_singletons == "FALSE"){
    detect_singletons = FALSE
}
if (detect_singletons == "true" || detect_singletons == "TRUE"){
    detect_singletons = TRUE
}

#Set DADA options
setDadaOpt(OMEGA_A = omegaa, OMEGA_P = omegap, OMEGA_C = omegac, DETECT_SINGLETONS = detect_singletons, BAND_SIZE = band_size)

#output_dir
path_results = Sys.getenv('output_dir')

### Denoise
cat(";; Working directory = ", workingDir)
cat(";; Performing DADA2 denoising ;; ")
#check for output dir and delete if needed
if (dir.exists(path_results)) {
    unlink(path_results, recursive=TRUE)
}
#create output dir
dir.create(path_results)

#filtered files path
qFiltered = readRDS(file.path(workingDir, "qFiltered.rds"))
sample_names = readRDS(file.path(workingDir, "sample_names.rds"))
cat("sample names = ", sample_names, "\n")
names(qFiltered) = sample_names

#Learn errors
set.seed(100)
errors = learnErrors(qFiltered,
            errorEstimationFunction = errorEstFun,
            BAND_SIZE = band_size,
            multithread = TRUE,
            verbose = TRUE,
            randomize = TRUE)
saveRDS(errors, file.path(path_results, "errors.rds"))

#Error rate figures
pdf(file.path(path_results, "Error_rates.pdf"))
    print( plotErrors(errors) )
dev.off()

# Infer sequence variants
denoised = vector("list", length(sample_names))
names(denoised) = sample_names
for(sam in sample_names) {
  cat("Processing:", sam, "\n")
  #Dereplicate
  dereplicated = derepFastq(qFiltered[[sam]],
                verbose = TRUE,
                qualityType = qualityType)
  #Denoise
  denoised[[sam]] = dada(dereplicated,
                err = errors,
                BAND_SIZE = band_size,
                multithread = TRUE,
                verbose = TRUE)
}
saveRDS(denoised, file.path(path_results, "denoised.rds"))


### WRITE PER-SAMPLE DENOISED and MERGED FASTA FILES
#make sequence table
cat(";; Writing per-sample denoised and merged fasta files")
ASV_tab = makeSequenceTable(denoised)
#write RDS object
saveRDS(ASV_tab, (file.path(path_results, "ASVs_table.denoised.rds")))

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

    #remove 0 seqs samples from qfilt statistics
    qfilt = readRDS(file.path(workingDir, "quality_filtered_read_count.rds"))
    row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
    qfilt = qfilt[row_sub, ]

seq_count <- cbind(qfilt, sapply(denoised, getN))
colnames(seq_count) <- c("input", "qualFiltered", "denoised")
write.csv(seq_count, 
            file.path(path_results, "seq_count_summary.csv"),
            row.names = TRUE,
            quote = FALSE)
cat(";; DONE")

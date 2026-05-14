#!/usr/bin/env Rscript

#DADA2 denoising single-end data.

#load dada2
library('dada2')

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

#load env variables
readType = Sys.getenv('readType')
fileFormat = Sys.getenv('fileFormat')
dataFormat = Sys.getenv('dataFormat')
workingDir = Sys.getenv('workingDir')

#load  variables
errorEstFun = Sys.getenv('errorEstFun')
pool = Sys.getenv('pool')
qualityType = Sys.getenv('qualityType')
nbases = as.numeric(Sys.getenv('nbases'))
randomize = Sys.getenv('randomize')

#setDadaOpt() settings
omegaa = as.numeric(Sys.getenv('OMEGA_A'))
omegap = as.numeric(Sys.getenv('OMEGA_P'))
omegac = as.numeric(Sys.getenv('OMEGA_C'))
detect_singletons = Sys.getenv('DETECT_SINGLETONS')
band_size = as.numeric(Sys.getenv('BAND_SIZE'))
homopoly_gap_penalty_raw = Sys.getenv('Homopoly_gap_penalty')

#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (homopoly_gap_penalty_raw %in% c("", "NULL", "null", "0")) {
    homopoly_gap_penalty = NULL
} else {
    homopoly_gap_penalty = as.numeric(homopoly_gap_penalty_raw)
}

cat(";; Settings:\n")
cat(";; errorEstimationFunction = ", errorEstFun, "\n")
cat(";; BAND_SIZE = ", band_size, "\n")
cat(";; OMEGA_A = ", omegaa, "\n")
cat(";; OMEGA_P = ", omegap, "\n")
cat(";; OMEGA_C = ", omegac, "\n")
cat(";; DETECT_SINGLETONS = ", detect_singletons, "\n")
cat(";; nbases = ", nbases, "\n")
cat(";; randomize = ", randomize, "\n")
cat(";; HOMOPOLYMER_GAP_PENALTY = ",
    if (is.null(homopoly_gap_penalty)) "NULL" else homopoly_gap_penalty, "\n")
cat(";; pool = ", pool, "\n")
cat(";; qualityType = ", qualityType, "\n")
cat(";; \n")

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
setDadaOpt(OMEGA_A = omegaa, OMEGA_P = omegap, OMEGA_C = omegac, DETECT_SINGLETONS = detect_singletons, BAND_SIZE = band_size, HOMOPOLYMER_GAP_PENALTY = homopoly_gap_penalty)

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
            nbases = nbases,
            BAND_SIZE = band_size,
            multithread = TRUE,
            verbose = TRUE,
            randomize = randomize)
saveRDS(errors, file.path(path_results, "errors.rds"))

#Error rate figures
pdf(file.path(path_results, "Error_rates.pdf"))
    print( plotErrors(errors) )
dev.off()

# Infer sequence variants
# pool=FALSE: dereplicate then denoise per sample (lowest peak RAM)
# pool=TRUE or "pseudo": dereplicate all samples into a named list, then one dada() on that list.
if (identical(pool, FALSE)) {
  cat(";; Denoising per sample with pool = ", pool, "\n", sep = "")
  denoised = vector("list", length(sample_names))
  names(denoised) = sample_names
  for (sam in sample_names) {
    cat(";; Dereplicating and denoising:", sam, "\n")
    dr = derepFastq(qFiltered[[sam]],
                  verbose = TRUE,
                  qualityType = qualityType)
    denoised[[sam]] = dada(dr,
                  err = errors,
                  pool = FALSE,
                  multithread = TRUE,
                  verbose = TRUE)
  }
} else {
  dereplicated = vector("list", length(sample_names))
  names(dereplicated) = sample_names
  for (sam in sample_names) {
    cat(";; Dereplicating:", sam, "\n")
    dereplicated[[sam]] = derepFastq(qFiltered[[sam]],
                  verbose = TRUE,
                  qualityType = qualityType)
  }
  cat(";; Denoising all samples with pool = ", pool, " (single dada() on list)\n", sep = "")
  denoised = dada(dereplicated,
                err = errors,
                pool = pool,
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

#!/usr/bin/env Rscript

#DADA2 denoising and merging paired-end data, for DADA2 full workflow.

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
minOverlap = as.numeric(Sys.getenv('minOverlap'))
maxMismatch = as.numeric(Sys.getenv('maxMismatch'))
trimOverhang = Sys.getenv('trimOverhang')
justConcatenate = Sys.getenv('justConcatenate')
pool = Sys.getenv('pool')
qualityType = Sys.getenv('qualityType')
errorEstFun = Sys.getenv('errorEstFun')
band_size = as.numeric(Sys.getenv('BAND_SIZE'))
#setDadaOpt() settings
omegaa = as.numeric(Sys.getenv('OMEGA_A')
omegap = as.numeric(Sys.getenv('OMEGA_P')
omegac= as.numeric(Sys.getenv('OMEGA_C')
detect_singletons = Sys.getenv('DETECT_SINGLETONS')

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

#output_dir
path_results = Sys.getenv('output_dir')

#Set DADA options
setDadaOpt(OMEGA_A = omegaa, OMEGA_P = omegap, OMEGA_C = omegac, DETECT_SINGLETONS = detect_singletons, BAND_SIZE = band_size)

### Denoise
if (pool != ""){
    cat("| Working directory = ", workingDir)
    cat("| Performing DADA2 denoising | ")
    cat("errorEstimationFunction = ", errorEstFun, "\n")
    cat("BAND_SIZE = ", band_size, "\n\n")  
    #check for output dir and delete if needed
    if (dir.exists(path_results)) {
        unlink(path_results, recursive=TRUE)
    }
    #create output dir
    dir.create(path_results)

    #filtered files path
    filtFs = readRDS(file.path(workingDir, "filtFs.rds"))
    filtRs = readRDS(file.path(workingDir, "filtRs.rds"))
    sample_names = readRDS(file.path(workingDir, "sample_names.rds"))

    #Learn the error rates
    errF = learnErrors(filtFs, errorEstimationFunction = errorEstFun, multithread = FALSE)
    errR = learnErrors(filtRs, errorEstimationFunction = errorEstFun, multithread = FALSE)

    #Error rate figures
    pdf(file.path(path_results, "Error_rates_R1.pdf"))
      print( plotErrors(errF) )
    dev.off()
    pdf(file.path(path_results, "Error_rates_R2.pdf"))
      print( plotErrors(errR) )
    dev.off()

    #dereplicate
    derepFs = derepFastq(filtFs, qualityType = qualityType)
    derepRs = derepFastq(filtRs, qualityType = qualityType)
    saveRDS(derepFs, (file.path(path_results, "derepFs.rds")))
    saveRDS(derepRs, (file.path(path_results, "derepRs.rds")))

    #denoise
    dadaFs = dada(derepFs, err = errF, pool = pool, multithread = FALSE)
    dadaRs = dada(derepRs, err = errR, pool = pool, multithread = FALSE)
    saveRDS(dadaFs, (file.path(path_results, "dadaFs.rds")))
    saveRDS(dadaRs, (file.path(path_results, "dadaRs.rds")))
}

### Merge denoised paired-end reads
if (pool == ""){
    cat("| Working directory = ", path_results)
    cat("| Merging data with mergePairs | ")
    #load denoised data
    dadaFs = readRDS(file.path(path_results, "dadaFs.rds"))
    dadaRs = readRDS(file.path(path_results, "dadaRs.rds"))
    derepFs = readRDS(file.path(path_results, "derepFs.rds"))
    derepRs = readRDS(file.path(path_results, "derepRs.rds"))

    #merge paired-end reads
    merge = mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                            maxMismatch = maxMismatch,
                            minOverlap = minOverlap,
                            justConcatenate = justConcatenate,
                            trimOverhang = trimOverhang)

    ### WRITE PER-SAMPLE DENOISED and MERGED FASTA FILES
    #make sequence table
    ASV_tab = makeSequenceTable(merge)
    rownames(ASV_tab) = gsub("_R1.*", "", rownames(ASV_tab))
    
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
    sample_names = readRDS(file.path(workingDir, "sample_names.rds"))
    qfilt = readRDS(file.path(workingDir, "quality_filtered.rds"))
    
        #remove 0 seqs samples from qfilt statistics
        row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
        qfilt = qfilt[row_sub, ]

    seq_count <- cbind(qfilt, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merge, getN))
    colnames(seq_count) <- c("input", "qualFiltered", "denoised_R1", "denoised_R2", "merged")
    rownames(seq_count) <- sample_names
    write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), row.names = TRUE, quote = FALSE)
    print("DONE")
}

#!/usr/bin/env Rscript

#DADA2 denoising and merging paired-end data with loessErrfun, for DADA2 full workflow.
# PacBioErrfun is not included here.

#load dada2
library("dada2")

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
errorEstFun = "loessErrfun"

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

#output_dir
path_results = Sys.getenv('output_dir')

#Set DADA options
setDadaOpt(OMEGA_A = omegaa, OMEGA_P = omegap, OMEGA_C = omegac, DETECT_SINGLETONS = detect_singletons, BAND_SIZE = band_size)

### Denoise
if (pool != ""){
    cat(";; Working directory = ", workingDir)
    cat(";; errorEstimationFunction = ", errorEstFun, "\n")
    cat(";; BAND_SIZE = ", band_size, "\n")  
    cat(";; Performing DADA2 denoising ;; ")
    #check for output dir and delete if needed
    if (dir.exists(path_results)) {
        unlink(path_results, recursive=TRUE)
    }
    #create output dir
    dir.create(path_results)

    #copy rds files to denoised_assembled.dada2 (for making seq_count_summary)
    file.copy(file.path(workingDir, "sample_names.rds"), path_results)
    file.copy(file.path(workingDir, "quality_filtered.rds"), path_results)

    #filtered files path
    filtFs = readRDS(file.path(workingDir, "filtFs.rds"))
    filtRs = readRDS(file.path(workingDir, "filtRs.rds"))
    sample_names = readRDS(file.path(workingDir, "sample_names.rds"))
    cat(";; sample names = ", sample_names, ";; ")
    names(filtFs) = sample_names
    names(filtRs) = sample_names


    # Denoise based on errorEstimationFunction
    if (errorEstFun == "loessErrfun"){
        set.seed(100)
        #Learn R1 error rates
        errF = learnErrors(filtFs, errorEstimationFunction = loessErrfun, nbases = 1e8, multithread = TRUE) # nolint
        saveRDS(errF, (file.path(path_results, "errF.rds")))
        #Learn R2 error rates
        errR = learnErrors(filtRs, errorEstimationFunction = loessErrfun, nbases = 1e8, multithread = TRUE, randomize = TRUE) # nolint
        saveRDS(errR, (file.path(path_results, "errR.rds")))
        #Error rate figures
        cat(";; ")
        pdf(file.path(path_results, "Error_rates_R1.pdf"))
          print( plotErrors(errF) )
        dev.off()
        pdf(file.path(path_results, "Error_rates_R2.pdf"))
          print( plotErrors(errR) )
        dev.off()

        # Sample inference and merger of paired-end reads
        mergers = vector("list", length(sample_names))
        names(mergers) = sample_names
        for(sample in sample_names) {
          cat(";; Processing:", sample, "\n")
          derepF = derepFastq(filtFs[[sample]])
          ddF = dada(derepF, err = errF, multithread = TRUE)
          derepR = derepFastq(filtRs[[sample]])
          ddR = dada(derepR, err = errR, multithread = TRUE)
          merger = mergePairs(ddF, derepF, ddR, derepR)
          mergers[[sample]] = merger
        }
        rm(derepF); rm(derepR)
        gc()
        saveRDS(mergers, (file.path(path_results, "mergers.rds")))
    }
}

### Merge denoised paired-end reads
if (pool == ""){
    cat(";; Working directory = ", workingDir)
    cat(";; Merging data with mergePairs ")
    #load denoised data
    mergers = readRDS(file.path(path_results, "mergers.rds"))

    ### WRITE PER-SAMPLE DENOISED and MERGED FASTA FILES
    #make sequence table
    cat(";; Writing per-sample denoised and merged fasta files ")
    ASV_tab = makeSequenceTable(mergers)
    #rownames(ASV_tab) = gsub("_R1.*", "", rownames(ASV_tab)) #no need when doing "names(mergers) = sample_names" above
    
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
    getN = function(x) sum(getUniques(x))
    sample_names = readRDS(file.path(workingDir, "sample_names.rds"))
    qfilt = readRDS(file.path(workingDir, "quality_filtered.rds"))

        #remove 0 seqs samples from qfilt statistics
        row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
        qfilt = qfilt[row_sub, ]

    #seq_count = cbind(qfilt, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merge, getN))
    #colnames(seq_count) = c("input", "qualFiltered", "denoised_R1", "denoised_R2", "merged")
    seq_count = cbind(qfilt, sapply(mergers, getN))
    colnames(seq_count) = c("input", "qualFiltered", "merged")
    rownames(seq_count) = sample_names
    write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), row.names = TRUE, quote = FALSE)

    unlink(c("sample_names.rds", "quality_filtered.rds")) #not needed here anymore, present in qualFiltered_out dir

    cat(";; DONE")
}

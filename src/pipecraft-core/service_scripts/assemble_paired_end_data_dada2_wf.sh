#!/bin/bash

# denoise and assemble paired-end data with DADA2 dada and mergePairs functions. For DADA2 full workflow.

##########################################################
###Third-party applications:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
##########################################################

#load env variables
readType=${readType}
dataFormat=${dataFormat}
workingDir=${workingDir}

#load variables
minOverlap=${minOverlap}
maxMismatch=${maxMismatch}
trimOverhang=${trimOverhang}
justConcatenate=${justConcatenate}
pool=${pool}
qualityType=${qualityType}

errorEstFun=${errorEstFun}
band_size=${BAND_SIZE}
omegaa=${OMEGA_A}
omegap=${OMEGA_P}
omegac=${OMEGA_C}
detect_singletons=${DETECT_SINGLETONS}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dirs
output_dir=$"/input/denoised_assembled.dada2"
export output_dir

### Check that at least 2 samples are provided
files=$(ls /input/qualFiltered_out | grep -c ".$fileFormat")
if (( $files < 4 )); then
    printf '%s\n' "ERROR]: please provide at least 2 samples for the ASVs workflow
>Quitting" >&2
    end_process
fi

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Process samples with dada2 removeBimeraDenovo function in R
printf "# Running DADA2 denoising and assembling \n"
Rlog=$(Rscript /scripts/submodules/dada2_denoise_assemble_wf.R 2>&1)
echo $Rlog >> $output_dir/denoise_assemble.log 
wait
#format R-log file
sed -i "s/;; /\n/g" $output_dir/denoise_assemble.log

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    $output_dir/denoise_assemble.log
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
if [[ ! -z $omegaa ]]; then
    printf "# Denoising and assembling of paired-end sequencing data was performed with dada2 (see 'Core commands' below for the used settings).

### NOTE: ### 
Input sequences must be made up only of A/C/G/T for denoising (i.e maxN must = 0 in quality filtering step). Otherwise DADA2 fails, and no output is generated.
#############

Files in 'denoised_assembled.dada2':
# *ASVs.fasta   = denoised and assembled sequences per sample in FASTA format (no fastq output). 'Size' denotes the abundance of the ASV sequence.  
# Error_rates_R1.pdf    = plots for estimated R1 error rates
# Error_rates_R2.pdf    = plots for estimated R2 error rates
# seq_count_summary.csv = summary of sequence counts per sample
# *.rds = R objects for dada2.

Core commands -> 
setDadaOpt(OMEGA_A = $OMEGA_A, OMEGA_P = $omegap, OMEGA_C = $omegac, DETECT_SINGLETONS = $detect_singletons, BAND_SIZE = $band_size)
learn errors: errF = learnErrors(fnFs, errorEstimationFunction = $errorEstFun)
              errR = learnErrors(fnRs, errorEstimationFunction = $errorEstFun)
dereplicate:  derepFs = derepFastq(fnFs, qualityType = $qualityType)
              derepRs = derepFastq(fnRs, qualityType = $qualityType)
denoise:      dadaFs = dada(derepFs, err = errF, pool = $pool)
              dadaRs = dada(derepRs, err = errR, pool = $pool)
              
    Total run time was $runtime sec for denoising.
    " > $output_dir/README.txt
fi
if [[ -z $omegaa ]]; then
    printf "
assemble:     mergePairs(dadaFs, derepFs, dadaRs, derepRs, maxMismatch = $maxMismatch, minOverlap = $minOverlap, justConcatenate = $justConcatenate, trimOverhang = $trimOverhang)

Total run time was $runtime sec for merging pairs.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
########################################################" >> $output_dir/README.txt
fi

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
if [[ -z $pool ]]; then
    echo "fileFormat=fasta"
else
    echo "fileFormat=$extension"
fi
if [[ -z $pool ]]; then
    echo "readType=single_end"
else
    echo "readType=paired_end"
fi
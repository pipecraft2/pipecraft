#!/bin/bash

# denoise single-end data with DADA2.

################################################
###Third-party applications:
#dada2 v1.28
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# DADA2 version: $dada2_version\n"

#load env variables
readType=${readType}
extension=${fileFormat}
dataFormat=${dataFormat}
workingDir=${workingDir}

#load variables
minOverlap=${minOverlap}
maxMismatch=${maxMismatch}
trimOverhang=${trimOverhang}
justConcatenate=${justConcatenate}
pool=${pool}
qualityType=${qualityType}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dirs
output_dir=$"/input/denoised.dada2"
export output_dir

### Check that at least 2 samples are provided
files=$(ls $workingDir | grep -c ".$fileFormat")
if (( $files < 2 )); then
    printf '%s\n' "ERROR]: please provide at least 2 samples for the ASVs workflow
>Quitting" >&2
    end_process
fi

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
### Process samples with dada2 removeBimeraDenovo function in R
printf "# Running DADA2 denoising \n"
Rlog=$(Rscript /scripts/submodules/dada2_SE_denoise.R 2>&1)
echo $Rlog >> $output_dir/denoise.log 
wait
#format R-log file
sed -i "s/;; /\n/g" $output_dir/denoise.log


#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    rm $output_dir/denoise.log
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file 
printf "# Single-end sequencing data was denoised with DADA2 (see 'Core command' below for the used settings).

### NOTE: ### 
Input sequences must be made up only of A/C/G/T for denoising (i.e maxN must = 0 in quality filtering step). Otherwise DADA2 fails, and no output is generated.
#############

Files in 'denoised_assembled.dada2':
# *.ASVs.fasta            = denoised and assembled ASVs per sample. 'Size' denotes the abundance of the ASV sequence.  
# Error_rates.pdf         = plots for estimated error rates
# seq_count_summary.csv   = summary of sequence and ASV counts per sample
# *.rds                   = R objects for dada2.

Core commands -> 
setDadaOpt(OMEGA_A = $OMEGA_A, OMEGA_P = $OMEGA_P, OMEGA_C = $OMEGA_C, DETECT_SINGLETONS = $DETECT_SINGLETONS, BAND_SIZE = $BAND_SIZE)
dereplicate:  dereplicated <- derepFastq(input, qualityType = $qualityType)
learn errors: errors = learnErrors(dereplicated, errorEstimationFunction = $errorEstFun, BAND_SIZE = $BAND_SIZE)
denoise:      denoised = dada(dereplicated, err = errors, BAND_SIZE = $BAND_SIZE)

Total run time was $runtime sec.
##############################################
###Third-party applications for this process:
#dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
##############################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"

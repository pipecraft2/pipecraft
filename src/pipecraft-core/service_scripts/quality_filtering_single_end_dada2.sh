#!/bin/bash

# Quality filter SINGLE-END sequencing data with dada2
# Input = single-end fastq files

##########################################################
###Third-party applications:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
##########################################################

#load env variables
readType=${readType}
extension=${fileFormat}
dataFormat=${dataFormat}
workingDir=${workingDir}

#load variables
samp_ID=${samp_ID}
maxEE=${maxEE}
maxN=${maxN}
truncQ=${truncQ}
truncLen_R1=${truncLen}
minLen=${minLen}
maxLen=${maxLen}
minQ=${minQ}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/qualFiltered_out"
export output_dir
#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env

### Process samples with dada2 filterAndTrim function in R
printf "# Running DADA2 filterAndTrim \n"
Rlog=$(Rscript /scripts/submodules/dada2_SE_filterAndTrim.R 2>&1)
echo $Rlog > $output_dir/qFilt.log 
wait
printf "\n DADA2 filterAndTrim completed \n"

#format R-log file
sed -i "s/;; /\n/g" $output_dir/qFilt.log 

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    rm $output_dir/qFilt.log 
fi

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering was performed using dada2 (see 'Core command' below for the used settings).

Files in 'qualFiltered_out':
    # *.$extension             = quality filtered sequences per sample.
    # seq_count_summary.csv    = summary of sequence counts per sample.
    # *.rds                    = R objects for dada2.

Core command -> 
filterAndTrim(inputR1, outputR1, maxN = $maxN, maxEE = $maxEE, truncQ = $truncQ, truncLen = $truncLen_R1, maxLen = $maxLen, minLen = $minLen, minQ = $minQ, rm.phix = TRUE)

Total run time was $runtime sec.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
########################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"


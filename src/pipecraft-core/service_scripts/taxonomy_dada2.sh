#!/bin/bash

#DADA2 RDP naive Bayesian classifier (function assignTaxonomy)
#Input = fasta file in the working directory and specified database file

################################################
###Third-party applications:
#dada2 v1.28
##############################################
# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# DADA2 version: $dada2_version\n"

#env variables
workingDir=${workingDir}
extension=$fileFormat && export fileFormat 
#load variables
minBoot=${minBoot}
tryRC=${tryRC}
dada2_database=${dada2_database}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/taxonomy_out.dada2"
export output_dir

#start time
start_time=$(date)
start=$(date +%s)

### Check if files with specified extension exist in the dir
first_file_check
### Check if single-end files are compressed (decompress and check)
check_gz_zip_SE

### Get input fasta
i=$"0"
for file in *.$extension; do
    input_fasta=$(echo $file)
    i=$((i + 1))
done
if [[ $i > 1 ]]; then
    if [[ -s $workingDir/ASVs_lenFilt.fasta ]] && [[ -s $workingDir/ASVs_collapsed.fasta ]]; then #if table filtering was performed by collapsing identical ASVs and by length
        input_fasta=$"/input/ASVs_lenFilt.fasta"
    else 
        printf '%s\n' "ERROR]: more than one representative sequence file ($extension file) in the working folder" >&2
        end_process
    fi
    printf "\n input fasta = $input_fasta \n"
else
    printf "\n input fasta = $input_fasta \n"
fi

#############################
### Start of the workflow ###
#############################
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env

###Run DADA2 classifier in R
printf "# Running DADA2 classifier \n"
Rlog=$(Rscript /scripts/submodules/dada2_classifier.R 2>&1)
echo $Rlog > $output_dir/dada2_classifier.log 
wait
printf "\n DADA2 classifier completed \n"

########################################
### CLEAN UP AND COMPILE README FILE ###
########################################
#Delete tempdir
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    if [[ -f $output_dir/dada2_classifier.log ]]; then
        rm -f $output_dir/dada2_classifier.log
    fi
fi

end=$(date +%s)
runtime=$((end-start))

###Make README.txt file
printf "# Taxonomy was assigned using DADA2 classifier (see 'Core command' below for the used settings).

Query    = $input_fasta
Database = $dada2_database

# taxonomy.csv = classifier results with bootstrap values.
  [if the above file does not exists, then check if the database formatting is appropriate for DADA2 classifier]

Core command -> 
assignTaxonomy($input_fasta, $dada2_database, minBoot = $minBoot, tryRC = $tryRC, outputBootstraps = TRUE)

Total run time was $runtime sec.

################################################
###Third-party applications:
#dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
##############################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

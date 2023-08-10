#!/bin/bash

# Quality filter PAIRED-END sequencing data with dada2
# Input = paired-end fastq files

##########################################################
###Third-party applications:
#dada2 v1.20
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
#seqkit v2.0.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
##########################################################

#load env variables
readType=${readType}
extension=${fileFormat}
dataFormat=${dataFormat}
workingDir=${workingDir}

#load variables
read_R1=${read_R1}
read_R2=${read_R2}
samp_ID=${samp_ID}
maxEE=${maxEE}
maxN=${maxN}
truncQ=${truncQ}
truncLen_R1=${truncLen}
truncLen_R2=${truncLen_R2}
minLen=${minLen}
maxLen=${maxLen}
minQ=${minQ}
matchIDs=${matchIDs}

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
### Prepare working env and check paired-end data
prepare_PE_env
### Check file formatting for FASTQ 
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]] || [[ $extension == "fastq.gz" ]] || [[ $extension == "fq.gz" ]]; then
    :
else
    printf '%s\n' "ERROR]: $file formatting not supported here!
Supported extensions: fastq, fq (and gz or zip compressed formats).
>Quitting" >&2
    end_process
fi

#Check identifiers
if [[ -z $read_R1 ]] || [[ -z $read_R2 ]] || [[ -z $samp_ID ]]; then
    printf '%s\n' "ERROR]: 'read R1/R2' or 'samp_ID' are not specified.
    >Quitting" >&2
    end_process
fi
read_R1_a=$(echo $read_R1 | sed -e 's/\\//') #if dot is the separator, then remove \ from the read identifier
while read file; do
    if [[ $file == *"$read_R1_a"* ]]; then
        :
    else
        printf '%s\n' "ERROR]: 'read R1/R2' identifiers are incorrectly specified.
        Check also 'samp ID' setting.
        >Quitting" >&2
        end_process
    fi
done < tempdir2/paired_end_files.txt

### Process samples with dada2 filterAndTrim function in R
printf "# Running DADA2 filterAndTrim \n"
Rlog=$(Rscript /scripts/submodules/dada2_PE_filterAndTrim.R 2>&1)
echo $Rlog > $output_dir/dada2_PE_filterAndTrim.log 
wait
printf "\n DADA2 filterAndTrim completed \n\n"


### Synchronizing R1 and R2 reads if $matchIDs == "true" - WORK WITH SEQKIT for matchIDs = TRUE, because sometimes DADA2 CANNOT automatically identify paired-end headers
if [[ $matchIDs == "true" ]] || [[ $matchIDs == "TRUE" ]]; then
    while read LINE; do
        #Read in R1 and R2 file names; without extension
        samp_name=$(basename $LINE | awk -F\\${samp_ID} '{print$1}')
        #If outputs are not empty, then synchronize R1 and R2
        if [[ -s $output_dir/$samp_name\_R1.$extension ]]; then
            if [[ -s $output_dir/$samp_name\_R2.$extension ]]; then
                printf "\nSynchronizing $samp_name R1 and R2 reads\n"
                cd $output_dir
                checkerror=$(seqkit pair -1 $samp_name\_R1.$extension -2 $samp_name\_R2.$extension 2>&1)
                check_app_error

                rm $samp_name\_R1.$extension
                rm $samp_name\_R2.$extension
                mv $samp_name\_R1.paired.$extension $samp_name\_R1.$extension
                mv $samp_name\_R2.paired.$extension $samp_name\_R2.$extension
                cd ..
            fi
        else
            printf "NOTE: all reads descarded from $samp_name\n"
        fi
    done < tempdir2/paired_end_files.txt
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ -d tempdir2 ]]; then
    rm -rf tempdir2
fi

### end pipe if no outputs were generated
outfile_check=$(ls $output_dir/*.$extension 2>/dev/null | wc -l)
if [[ $outfile_check != 0 ]]; then 
    qFilt_compleated=$"TRUE"
    export qFilt_compleated
else 
    printf '%s\n' "ERROR]: no output files generated after quality filtering ($output_dir). Adjust settings.
    >Quitting" >&2
    end_process
fi

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering with dada2.

Files in 'qualFiltered_out':
    # *.$extension             = quality filtered sequences per sample.
    # seq_count_summary.txt    = summary of sequence counts per sample.
    # *.rds                    = R objects for dada2.

Core command -> 
filterAndTrim(inputR1, outputR1, inputR2, outputR2, maxN = $maxN, maxEE = c($maxEE, $maxEE), truncQ = $truncQ, truncLen= c($truncLen_R1, $truncLen_R2), maxLen = $maxLen, minLen = $minLen, minQ=$minQ, rm.phix = TRUE, multithread = TRUE)

Total run time was $runtime sec.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.20
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
#seqkit v2.0.0 for synchronizing R1 and R2 after filtering (when matchIDs = TRUE)
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
########################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "dataFormat=$dataFormat"
echo "readType=paired_end"


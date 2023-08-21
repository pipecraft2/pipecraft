#!/bin/bash

# ASSEMBLE PAIRED-END data with dada2
# Input = paired-end fastq files.

##########################################################
###Third-party applications:
#dada2 v1.26
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
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
read_R1=${read_R1}
read_R2=${read_R2}
samp_ID=${samp_ID}
minOverlap=${minOverlap}
maxMismatch=${maxMismatch}
trimOverhang=${trimOverhang}
justConcatenate=${justConcatenate}
pool=${pool}
qualityType=${qualityType}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/denoised_assembled.dada2"

### Check that at least 2 samples are provided
files=$(ls $workingDir | grep -c "$extension")
if [[ $files < 4 ]]; then
    printf '%s\n' "ERROR]: please provide at least 2 samples for the ASVs workflow
>Quitting" >&2
    end_process
fi

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env

#Check file formatting for FASTQ
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
while read file; do
    if [[ $file == *"$read_R1"* ]]; then
        :
    else
        printf '%s\n' "ERROR]: 'read R1/R2' identifiers are incorrectly specified.
        Check also 'samp ID' setting.
        >Quitting" >&2
        end_process
    fi
done < tempdir2/paired_end_files.txt

#Check Ns in the fastq files (not allowed for DADA2 denoising)
while read file; do
    find_Ns=$(seqkit grep --quiet --by-seq --pattern 'N' $file | wc -l)
    if [[ $find_Ns == 0 ]]; then
        :
    else
        printf '%s\n' "ERROR]: sequences must be made up only of A/C/G/T. Supply quality filtered fastq files.
        >Quitting" >&2
        end_process
    fi
done < tempdir2/files_in_folder.txt


### Process samples with dada2 filterAndTrim function in R
printf "# Running DADA2 mergePairs \n"
Rlog=$(Rscript /scripts/submodules/dada2_mergePairs.R 2>&1)
echo $Rlog > $output_dir/denoise_assemble.log 
wait
printf "\n DADA2 mergePairs completed \n"

# Rereplicate sequences per sample
for file in $output_dir/*.fasta; do
    samp_name=$(basename $file | awk -F "$samp_ID" '{print $1}')
    echo $samp_name
    vsearch --rereplicate $file --fasta_width 0 --output $output_dir/$samp_name.fasta -relabel $samp_name.
    rm $file
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ -d tempdir2 ]]; then
    rm -rf tempdir2
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Denoising and assembling of PAIRED-END sequencing data with dada2.

### NOTE: ### 
Input sequences must be made up only of A/C/G/T for denoising (i.e maxN must = 0 in quality filtering step). Otherwise DADA2 fails, and no output is generated.
#############

Files in 'denoised_assembled.dada2':
# *ASVs.fasta   = denoised and assembled sequences per sample in FASTA format (no fastq output). 'Size' denotes the abundance of the ASV sequence.  
# Error_rates_R1.pdf    = plots for estimated R1 error rates
# Error_rates_R2.pdf    = plots for estimated R2 error rates
# seq_count_summary.txt = summary of sequence counts per sample
# *.rds = R objects for dada2.

Core commands -> 
learn errors: err = learnErrors(input)
dereplicate:  derep = derepFastq(input, qualityType = $qualityType)
denoise:      dadaFs = dada(input, err = err, pool = $pool)
assemble:     mergePairs(inputR1, dereplicatedR1, inputR2, dereplicatedR2, maxMismatch = $maxMismatch, minOverlap = $minOverlap, justConcatenate = $justConcatenate, trimOverhang = $trimOverhang)

Total run time was $runtime sec.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.26
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
########################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=single_end"


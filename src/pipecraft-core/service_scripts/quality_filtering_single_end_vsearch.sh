#!/bin/bash

# Quality filter SINGLE-END sequencing data with vsearch
# Input = single-end fastq files

##########################################################
###Third-party applications:
#vsearch v2.23.0
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz v2.4
##########################################################

#load variables
extension=$fileFormat && export fileFormat 
maxee=$"--fastq_maxee ${maxee}"
maxns=$"--fastq_maxns ${maxNs}"
minlen=$"--fastq_minlen ${min_length}"
cores=$"--threads ${cores}"
qmax=$"--fastq_qmax ${qmax}"
qmin=$"--fastq_qmin ${qmin}"
trunc_length=$trunc_length
max_length=$max_length
maxee_rate=$maxee_rate
truncqual=${truncqual}
truncee=${truncee}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/qualFiltered_out"

#additional options, if selection != undefined
if [[ $max_length == null ]] || [[ -z $max_length ]]; then
    max_length=$""
else
    max_length=$"--fastq_maxlen $max_length"
fi
if [[ $maxee_rate == null ]] || [[ -z $maxee_rate ]]; then
    maxee_rate=$""
else
    maxee_rate=$"--fastq_maxee_rate $maxee_rate"
fi
if [[ $trunc_length == null ]] || [[ -z $trunc_length ]]; then
    trunc_length=$""
else
    trunc_length=$"--fastq_trunclen $trunc_length"
fi
if [[ $truncqual == null ]] || [[ -z $truncqual ]]; then
    truncqual=$""
else
    truncqual=$"--fastq_truncqual $truncqual"
fi
if [[ $truncee == null ]] || [[ -z $truncee ]]; then
    truncee=$""
else
    truncee=$"--fastq_truncee $truncee"
fi

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
### Process samples
for file in *.$extension; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$extension//")
    ## Preparing files for the process
    printf "\n____________________________________\n"
    printf "Processing $input ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        
    check_gz_zip_SE
    ### Check input formats (fastq supported)
    check_extension_fastq

    ###############################
    ### Start quality filtering ###
    ###############################
    checkerror=$(vsearch --fastq_filter \
    $input.$extension \
    $maxee \
    $maxns \
    $trunc_length \
    $minlen \
    $cores \
    $qmax \
    $qmin \
    $max_length \
    $maxee_rate \
    $truncqual \
    $truncee \
    --fastqout $output_dir/$input.$extension 2>&1)
    check_app_error
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering was performed using vsearch (see 'Core command' below for the used settings).

Files in 'qualFiltered_out':
# *.$extension           = quality filtered sequences in FASTQ format.
# seq_count_summary.txt     = summary of sequence counts per sample.

Core command -> 
vsearch --fastq_filter input_file $maxee $maxns $trunc_length $minlen $cores $qmax $qmin $max_length $maxee_rate $truncqual $truncee --fastqout $output_dir/output_file.fastq

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.23.0 for quality filtering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
##########################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

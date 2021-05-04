#!/bin/bash

#Input = single-end fastq files.

# Quality filter SINGLE-END sequencing data with trimmomatic

##########################################################
###Third-party applications:
#trimmomatic
    #citation: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu1
    #Distributed under the GNU GENERAL PUBLIC LICENE
    #https://github.com/usadellab/Trimmomatic
#seqkit
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz
##########################################################

###############################
###############################
#These variables are for testing (DELETE when implementing to PipeCraft)
extension=$"fastq"
#mandatory options
window_size=$"5"
required_qual=$"27"
min_length=$"32"
#additional options
threads=$"4"
phred="33"
leading_qual_threshold=$"11" #or 'undefined', if selection is not active
trailing_qual_threshold=$"11" #or 'undefined', if selection is not active
###############################
###############################

#############################
### Start of the workflow ###
#############################
#additional options, if selection != undefined
if [[ $leading_qual_threshold == "undefined" ]]; then
    :
else
    LEADING=$"LEADING:$leading_qual_threshold"
fi
if [[ $trailing_qual_threshold == "undefined" ]]; then
    :
else
    TRAILING=$"TRAILING:$trailing_qual_threshold"
fi

start=$(date +%s)
# Source for functions
source /scripts/framework.functions.sh
#output dir
output_dir=$"qualFiltered_out"
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_SE_env
### Process samples
for file in *.$extension; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$extension//")
    ## Preparing files for the process
    printf "\n____________________________________\n"
    printf "Processing $input ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        #$extension will be $newextension
    check_gz_zip_SE
    ### Check input formats (fastq supported)
    check_extension_fastq

    ###############################
    ### Start quality filtering ###
    ###############################
    checkerror=$(trimmomatic SE \
    $input.$newextension \
    $output_dir/$input.qualFilt.$newextension \
    -phred$phred \
    $LEADING \
    $TRAILING \
    SLIDINGWINDOW:$window_size:$required_qual \
    MINLEN:$min_length \
    -threads $threads 2>&1)
    check_app_error

    #Convert output fastq files to FASTA
    mkdir -p $output_dir/qualFilt_FASTA
    checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$input.qualFilt.$newextension -o $output_dir/FASTA/$input.qualFilt.fasta 2>&1)
    check_app_error
done



#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
#file identifier string after the process
outfile_addition=$"qualFilt"
clean_and_make_stats

#Make README.txt file
printf "Files in /$output_dir directory represent quality filtered sequences in FASTQ format according to the selected options.
Files in $output_dir/FASTA directory represent quality filtered sequences in FASTA format.
If the quality of the data is sufficent after this step (check with FastQC module), then
you may proceed with FASTA files only.\n" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"

end=$(date +%s)
runtime=$((end-start))
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=/$output_dir"
echo "fileFormat=$newextension"
echo "dataFormat=$dataFormat"
echo "readType=single-end"
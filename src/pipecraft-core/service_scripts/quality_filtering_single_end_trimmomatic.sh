#!/bin/bash

# Quality filter SINGLE-END sequencing data with trimmomatic
#Input = single-end fastq files.

################################################
###Third-party applications:
#trimmomatic v0.39
    #citation: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu1
    #Distributed under the GNU GENERAL PUBLIC LICENE
    #https://github.com/usadellab/Trimmomatic
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz v2.4
################################################

#load variables
window_size=${window_size}
required_qual=${required_quality}
min_length=${min_length}
threads=${cores}
phred=${phred}
leading_qual_threshold=${leading_qual_threshold}
trailing_qual_threshold=${trailing_qual_threshold}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/qualFiltered_out"

#additional options, if selection != undefined
if [[ $leading_qual_threshold == null ]] || [[ -z $leading_qual_threshold ]]; then
    LEADING=$""
else
    LEADING=$"LEADING:$leading_qual_threshold"
fi
if [[ $trailing_qual_threshold == null ]] || [[ -z $trailing_qual_threshold ]]; then
    TRAILING=$""
else
    TRAILING=$"TRAILING:$trailing_qual_threshold"
fi

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
### Process samples
for file in *.$fileFormat; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
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
    checkerror=$(java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE \
    $input.$extension \
    $output_dir/$input.$extension \
    -phred$phred \
    $LEADING \
    $TRAILING \
    SLIDINGWINDOW:$window_size:$required_qual \
    MINLEN:$min_length \
    -threads $threads 2>&1)
    check_app_error

    #Convert output fastq files to FASTA
    mkdir -p $output_dir/FASTA
    checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$input.$extension -o $output_dir/FASTA/$input.fasta 2>&1)
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
printf "# Quality filtering was performed using trimmomatic (see 'Core commands' below for the used settings).

Files in 'qualFiltered_out':
# *.$extension              = quality filtered sequences in FASTQ format.
# seq_count_summary.txt     = summary of sequence counts per sample.
Files in 'qualFiltered_out/FASTA':
# *.fasta                   = quality filtered sequences in FASTA format.

Core commands ->
quality filtering: trimmomatic-0.39.jar SE input_file output_file -phred$phred $LEADING $TRAILING SLIDINGWINDOW:$window_size:$required_qual MINLEN:$min_length -threads $threads
convert output fastq files to FASTA: seqkit fq2fa -t dna --line-width 0 input_file -o FASTA/output_file.fasta

Total run time was $runtime sec.

##############################################
###Third-party applications for this process:
#trimmomatic v0.39 for quality filtering
    #citation: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu1
    #https://github.com/usadellab/Trimmomatic
#seqkit v2.3.0 for converting filtered fastq to fasta 
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
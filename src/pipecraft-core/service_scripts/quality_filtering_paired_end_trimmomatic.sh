#!/bin/bash

#Input = paired-end fastq files.

# Quality filter PAIRED-END sequencing data with trimmomatic

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
### Prepare working env and check paired-end data
prepare_PE_env
### Process samples
while read LINE; do
    #Read in R1 and R2 file names; without extension
    inputR1=$(echo $LINE | sed -e "s/.$fileFormat//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
    ## Preparing files for the process
    printf "\n____________________________________\n"
    printf "Processing $inputR1 and $inputR2 ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_gz_zip_PE
    ### Check input formats (fastq supported)
    check_extension_fastq

    ###############################
    ### Start quality filtering ###
    ###############################
    #make dir for discarded seqs
    mkdir -p $output_dir/discarded

    checkerror=$(java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE \
    $inputR1.$extension $inputR2.$extension \
    $output_dir/$inputR1.$extension $output_dir/discarded/$inputR1.discarded.$extension \
    $output_dir/$inputR2.$extension $output_dir/discarded/$inputR2.discarded.$extension \
    $LEADING \
    $TRAILING \
    -phred$phred \
    SLIDINGWINDOW:$window_size:$required_qual \
    MINLEN:$min_length \
    -threads $threads 2>&1)
    check_app_error

    #Convert output fastq files to FASTA
    mkdir -p $output_dir/FASTA
    checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$inputR1.$extension -o $output_dir/FASTA/$inputR1.fasta 2>&1)
    check_app_error
    checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$inputR2.$extension -o $output_dir/FASTA/$inputR2.fasta 2>&1)
    check_app_error
done < tempdir2/paired_end_files.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for discarded seqs
printf "Files in /discarded folder represent sequences that did not pass quality filtering.\n
If no files in this folder, then all sequences were passed to files in $output_dir directory" > $output_dir/untrimmed/README.txt

#Make README.txt file
printf "# Quality filtering was performed using trimmomatic (see 'Core commands' below for the used settings).

Files in 'qualFiltered_out':
# *.$extension              = quality filtered sequences in FASTQ format.
# seq_count_summary.txt     = summary of sequence counts per sample.

Files in 'qualFiltered_out/FASTA':
# *.fasta                   = quality filtered sequences in FASTA format.
Files in 'qualFiltered_out/discarded':
# *.discarded.$extension = discarded sequences.

Core commands -> 
quality filtering: trimmomatic-0.39.jar PE inputR1 inputR2 outputR1 discarded/outputR1.discarded outputR2 discarded/outputR2.discarded $LEADING $TRAILING -phred$phred SLIDINGWINDOW:$window_size:$required_qual MINLEN:$min_length -threads $threads
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
echo "readType=paired_end"

#!/bin/bash

# Quality filter PAIRED-END sequencing data with vsearch
# Input = paired-end fastq files

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
########################################################

#load variables
extension=$fileFormat
maxee=$"--fastq_maxee ${maxee}"
maxns=$"--fastq_maxns ${maxNs}"
minlen=$"--fastq_minlen ${min_length}"
cores=$"--threads ${cores}"
qmax=$"--fastq_qmax ${qmax}"
qmin=$"--fastq_qmin ${qmin}"
trunc_length=${trunc_length}
max_length=${max_length}
maxee_rate=${maxee_rate}
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
### Prepare working env and check paired-end data
prepare_PE_env
### Process samples
while read LINE; do
    #Read in R1 and R2 file names; without extension
    inputR1=$(echo $LINE | sed -e "s/.$extension//")
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
    mkdir -p tempdir
 
    printf "vsearch --fastq_filter \
    $inputR1.$extension \
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
    --fastqout tempdir/$inputR1.$extension"

    #R1
    checkerror=$(vsearch --fastq_filter \
    $inputR1.$extension \
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
    --fastqout tempdir/$inputR1.$extension 2>&1)
    check_app_error

    #R2
    checkerror=$(vsearch --fastq_filter \
    $inputR2.$extension \
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
    --fastqout tempdir/$inputR2.$extension 2>&1)
    check_app_error

    #If outputs are not empty, then synchronize R1 and R2
    if [ -s tempdir/$inputR1.$extension ]; then
        if [ -s tempdir/$inputR2.$extension ]; then
            printf "\nSynchronizing R1 and R2 reads (matching order for paired-end reads merging)\n"
            cd tempdir
            checkerror=$(seqkit pair -1 $inputR1.$extension -2 $inputR2.$extension 2>&1)
            check_app_error

            rm $inputR1.$extension
            rm $inputR2.$extension
            mv $inputR1.paired.$extension $inputR1.$extension
            mv $inputR2.paired.$extension $inputR2.$extension
            cd ..

            #Move final files to $output_dir
            mv tempdir/$inputR1.$extension $output_dir/$inputR1.$extension
            mv tempdir/$inputR2.$extension $output_dir/$inputR2.$extension

            # #Convert output fastq files to FASTA
            # mkdir -p $output_dir/FASTA
            # checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$inputR1.$extension -o $output_dir/FASTA/$inputR1.fasta 2>&1)
            # check_app_error
            # checkerror=$(seqkit fq2fa -t dna --line-width 0 $output_dir/$inputR2.$extension -o $output_dir/FASTA/$inputR2.fasta 2>&1)
            # check_app_error
        else
            printf '%s\n' "WARNING]: $inputR2 has 0 seqs after filtering (no output for that sample)"
            rm tempdir/$inputR1.$extension
            rm tempdir/$inputR2.$extension
        fi
    else
        printf '%s\n' "WARNING]: $inputR1 has 0 seqs after filtering (no output for that sample)"
        rm tempdir/$inputR1.$extension
        rm tempdir/$inputR2.$extension
    fi
done < tempdir2/paired_end_files.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering with vsearch.

Files in 'qualFiltered_out':
# *.$extension           = quality filtered sequences in FASTQ format.
# seq_count_summary.txt     = summary of sequence counts per sample.

Core commands -> 
quality filtering: vsearch --fastq_filter input_file $maxee $maxns $trunc_length $minlen $cores $qmax $qmin $max_length $maxee_rate $truncqual $truncee --fastqout output_file
Synchronizing R1 and R2 reads: seqkit pair -1 inputR1 -2 inputR2

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.23.0 for quality filtering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
#seqkit v2.3.0 for synchronizing R1 and R2 after filtering
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

echo "readType=paired_end"

#!/bin/bash

#Input = paired-end fastq or paired-end fasta files.

# REMOVE PRIMERS from paired-end reads (e.g. Illumina reads)
#Degenerate primers are allowed using IUPAC codes.

##########################################################
###Third-party applications:
#cutadapt
    #citation: Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), 10-12.
    #Distributed under MIT License"
    #https://cutadapt.readthedocs.io/en/stable/#
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
# echo $fileFormat
# echo $mismatches
# echo $min_seq_length
# echo $overlap
# echo $cores
# echo $no_indels
# echo $discard_untrimmed
# echo $seqs_to_keep
# echo $forward_primers
# echo $reverse_primers
# echo "-e $mismatches"
extension=$"fastq"
mismatches=$"-e 2"
min_length=$"--minimum-length 19"
overlap=$"--overlap 15"
cores=$"--cores 0"
no_indels=$"TRUE"
discard_untrimmed=$"TRUE"
seqs_to_keep=$"keep_all" #keep_all/keep_only_linked


# fwd_tempprimer=$forward_primers
# rev_tempprimer=$reverse_primers
fwd_tempprimer=$"ACCTGCGGARGGATCA"
rev_tempprimer=$"GAGATCCRTTGYTRAAAGTT"
# fwd_tempprimer=$"ACCTGCTAGGCTAGATGC,GCTAGCTAGCTAGCTGATGC,ATCGATGCTAGCTAGCTAGCTGA"
# rev_tempprimer=$"GGGATCCATCGATTTAAC,GCTAGCTAGCTAGCTAGCTAGC"

###############################
###############################

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
# Source for functions
source /scripts/framework.functions.sh
#output dir
output_dir=$"/input/primersCut_out"
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env

### Prepare primers
# Primer arrays
fwd_primer_array=$(echo $fwd_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
rev_primer_array=$(echo $rev_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
# Forward primer(s) to fasta file

i=1
for primer in $fwd_primer_array; do
    echo ">fwd_primer$i" >> tempdir2/fwd_primer.fasta
    echo $primer >> tempdir2/fwd_primer.fasta
    ((i=i+1))
done
# Reverse primer(s) to fasta file
i=1
for primer in $rev_primer_array; do
    echo ">rev_primer$i" >> tempdir2/rev_primer.fasta
    echo $primer >> tempdir2/rev_primer.fasta
    ((i=i+1))
done
# Reverse complement FWD and REV primers
checkerror=$( seqkit seq --quiet -t dna -r -p tempdir2/fwd_primer.fasta >> tempdir2/fwd_primer_RC.fasta 2>&1)
check_app_error
checkerror=$( seqkit seq --quiet -t dna -r -p tempdir2/rev_primer.fasta >> tempdir2/rev_primer_RC.fasta 2>&1)
check_app_error
# Make linked primers files
i=1
while read LINE; do
    fwd_primer=$(echo $LINE | grep -v ">")
    if [ -z "$fwd_primer" ]; then
        :
    else
        while read LINE; do
            rev_primer_RC=$(echo $LINE | grep -v ">")
            if [ -z "$rev_primer_RC" ]; then
                :
            else
                echo ">primer$i" >> tempdir2/liked_fwd_revRC.fasta
                echo "$fwd_primer...$rev_primer_RC" >> tempdir2/liked_fwd_revRC.fasta
                ((i=i+1))
            fi
        done < tempdir2/rev_primer_RC.fasta
    fi
done < tempdir2/fwd_primer.fasta
i=1
while read LINE; do
    rev_primer=$(echo $LINE | grep -v ">")
    if [ -z "$rev_primer" ]; then
        :
    else
        while read LINE; do
            fwd_primer_RC=$(echo $LINE | grep -v ">")
            if [ -z "$fwd_primer_RC" ]; then
                :
            else
                echo ">primer$i" >> tempdir2/liked_rev_fwdRC.fasta
                echo "$rev_primer...$fwd_primer_RC" >> tempdir2/liked_rev_fwdRC.fasta
                ((i=i+1))
            fi
        done < tempdir2/fwd_primer_RC.fasta
    fi
done < tempdir2/rev_primer.fasta

#############################
### Start of the workflow ###
#############################
if [[ $no_indels == "TRUE" ]]; then
    indels=$"--no-indels"
fi
### Read through each file in paired_end_files.txt
while read LINE; do
    #Write file name without extension
    inputR1=$(echo $LINE | sed -e "s/.$extension//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
    ## Preparing files
    printf "\n____________________________________\n"
    printf "Checking $inputR1 and $inputR2  ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        #$extension will be $newextension
    check_gz_zip_PE
    ### Check input formats (only fastq/fasta supported)
    check_extension_fastx

    ##############################
    ### Start clipping primers ###
    ##############################
    printf "\n# Clipping primers from $inputR1 and $inputR2 (paired-end reads). \n"
    printf " forward primer(s): $fwd_tempprimer\n"
    printf " reverse primer(s): $rev_tempprimer\n"

    #If discard_untrimmed = TRUE, then assigns outputs and make outdir
    if [[ $discard_untrimmed == "TRUE" ]]; then
        mkdir -p $output_dir/untrimmed
        untrimmed_output=$"--untrimmed-output $output_dir/untrimmed/$inputR1.untrimmed.$newextension"
        untrimmed_paired_output=$"--untrimmed-paired-output $output_dir/untrimmed/$inputR2.untrimmed.$newextension"
    fi

    ### Clip primers with cutadapt
    if [[ $seqs_to_keep == "keep_all" ]]; then
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=both \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/liked_rev_fwdRC.fasta \
        -g file:tempdir2/fwd_primer.fasta \
        -a file:tempdir2/rev_primer_RC.fasta \
        -g file:tempdir2/rev_primer.fasta \
        -a file:tempdir2/fwd_primer_RC.fasta \
        -G file:tempdir2/liked_fwd_revRC.fasta \
        -G file:tempdir2/liked_rev_fwdRC.fasta \
        -G file:tempdir2/rev_primer.fasta \
        -A file:tempdir2/fwd_primer_RC.fasta \
        -G file:tempdir2/fwd_primer.fasta \
        -A file:tempdir2/rev_primer_RC.fasta \
        -o $output_dir/$inputR1.primersCut.$newextension \
        -p $output_dir/$inputR2.primersCut.$newextension \
        $inputR1.$newextension $inputR2.$newextension 2>&1)
        check_app_error

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=both \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/liked_rev_fwdRC.fasta \
        -G file:tempdir2/liked_fwd_revRC.fasta \
        -G file:tempdir2/liked_rev_fwdRC.fasta \
        -o $output_dir/$inputR1.primersCut.$newextension \
        -p $output_dir/$inputR2.primersCut.$newextension \
        $inputR1.$newextension $inputR2.$newextension 2>&1)
        check_app_error
    fi
done < tempdir2/paired_end_files.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
#file identifier string after the process
outfile_addition=$"primersCut"
clean_and_make_stats

#Make README.txt file for untrimmed seqs
if [[ $discard_untrimmed == "TRUE" ]]; then
    printf "Files in /untrimmed folder represent sequences that did not contain specified primer strings.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
If no files in this folder, then all sequences were passed to files in $output_dir directory" > $output_dir/untrimmed/README.txt
fi

#Make README.txt file for PrimerClipped reads
printf "Files in /$output_dir folder represent sequences from where the PCR primers were recognized and clipped.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
Output R1 and R2 reads are synchronized for merging paired-end data. 
If no outputs were generated into /$output_dir, check your specified primer stings and adjust settings" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"

end=$(date +%s)
runtime=$((end-start))
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$newextension"
echo "dataFormat=$dataFormat"
echo "readType=paired-end"
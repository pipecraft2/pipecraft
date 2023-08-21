#!/bin/bash

# REMOVE PRIMERS from paired-end reads
# Degenerate primers are allowed using IUPAC codes. Reverse complementary strings will be also searched.
# Input = paired-end fastq or paired-end fasta files. If using fasta, then cores must = 1

##########################################################
###Third-party applications:
#cutadapt v4.4
    #citation: Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), 10-12.
    #Distributed under MIT License
    #https://cutadapt.readthedocs.io/en/stable/#
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz v2.4
##########################################################

#load variables
extension=$fileFormat && export fileFormat 
mismatches=$"-e ${mismatches}"
min_length=$"--minimum-length ${min_seq_length}"
overlap=$"--overlap ${min_overlap}"
cores=$"--cores ${cores}"
no_indels=$no_indels
discard_untrimmed=$"TRUE" #currently only fixed to TRUE 
seqs_to_keep=$seqs_to_keep
pair_filter=$pair_filter

fwd_tempprimer=$forward_primers
rev_tempprimer=$reverse_primers

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/primersCut_out"

#############################
### Start of the workflow ###
#############################
if [[ $no_indels == "TRUE" ]]; then
    indels=$"--no-indels"
fi
start=$(date +%s)
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
checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/fwd_primer.fasta >> tempdir2/fwd_primer_RC.fasta 2>&1)
check_app_error
checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/rev_primer.fasta >> tempdir2/rev_primer_RC.fasta 2>&1)
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


### Read through each file in paired_end_files.txt
while read LINE; do
    #Write file name without extension
    inputR1=$(echo $LINE | sed -e "s/.$extension//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
    ## Preparing files
    ### Check input formats (fastq/fasta/gz supported)
    check_extension_fastxgz

    ##############################
    ### Start clipping primers ###
    ##############################
    printf "\n# Clipping primers from $inputR1 and $inputR2 (paired-end reads). \n"
    printf " forward primer(s): $fwd_tempprimer\n"
    printf " reverse primer(s): $rev_tempprimer\n"

    #If discard_untrimmed = TRUE, then assigns outputs and make outdir
    if [[ $discard_untrimmed == "TRUE" ]]; then
        mkdir -p $output_dir/untrimmed
        untrimmed_output=$"--untrimmed-output $output_dir/untrimmed/$inputR1.untrimmed.$extension"
        untrimmed_paired_output=$"--untrimmed-paired-output $output_dir/untrimmed/$inputR2.untrimmed.$extension"
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
        --pair-filter=$pair_filter \
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
        -o $output_dir/$inputR1.round1.$extension \
        -p $output_dir/$inputR2.round1.$extension \
        $inputR1.$extension $inputR2.$extension 2>&1)
        check_app_error

        #round2; clipping if present, but no discarding
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        --pair-filter=$pair_filter \
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
        -o $output_dir/$inputR1.$extension \
        -p $output_dir/$inputR2.$extension \
        $output_dir/$inputR1.round1.$extension $output_dir/$inputR2.round1.$extension 2>&1)
        check_app_error
        rm $output_dir/$inputR1.round1.$extension $output_dir/$inputR2.round1.$extension

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=$pair_filter \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/liked_rev_fwdRC.fasta \
        -G file:tempdir2/liked_fwd_revRC.fasta \
        -G file:tempdir2/liked_rev_fwdRC.fasta \
        -o $output_dir/$inputR1.round1.$extension \
        -p $output_dir/$inputR2.round1.$extension \
        $inputR1.$extension $inputR2.$extension 2>&1)
        check_app_error

        #additional check of the primer presence
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        -g file:tempdir2/fwd_primer.fasta \
        -g file:tempdir2/fwd_primer_RC.fasta \
        -g file:tempdir2/rev_primer.fasta \
        -g file:tempdir2/rev_primer_RC.fasta \
        -G file:tempdir2/fwd_primer.fasta \
        -G file:tempdir2/fwd_primer_RC.fasta \
        -G file:tempdir2/rev_primer.fasta \
        -G file:tempdir2/rev_primer_RC.fasta \
        -o $output_dir/$inputR1.$extension \
        -p $output_dir/$inputR2.$extension \
        $output_dir/$inputR1.round1.$extension $output_dir/$inputR2.round1.$extension 2>&1)
        check_app_error
        rm $output_dir/$inputR1.round1.$extension
        rm $output_dir/$inputR2.round1.$extension
    fi
done < tempdir2/paired_end_files.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for untrimmed seqs
if [[ $discard_untrimmed == "TRUE" ]]; then
    printf "Files in /untrimmed folder represent sequences that did not contain specified primer strings.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
If no files in this folder, then all sequences were passed to files in $output_dir directory" > $output_dir/untrimmed/README.txt
fi

#Make README.txt file for PrimerClipped reads
printf "Files in 'primersCut_out' folder represent sequences from where the PCR primers were recognized and clipped.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
Output R1 and R2 reads are synchronized for merging paired-end data. 
If no outputs were generated into /$output_dir, check your specified primer strings and adjust settings.
\nSummary of sequence counts in 'seq_count_summary.txt'\n
\n\nTotal run time was $runtime sec.\n\n\n
##########################################################
###Third-party applications used for this process [PLEASE CITE]:
#cutadapt v4.4 for cutting the primers
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0 for generating reverse complementary primer strings
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in 'seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=paired_end"
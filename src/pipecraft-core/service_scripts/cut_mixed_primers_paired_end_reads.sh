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
extension=${fileFormat}  # KEEP THIS (removed in some other scripts)
mismatches=$"-e ${mismatches}"
min_length=$"--minimum-length 32"   # minimum len of the output sequence. FIXED to 32 (in order to avoid 0 len seqs) because cutadapt --minimum-length does not behave as expected
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
# if keep_only_linked, then linked primers are REQUIRED (default = optional)
if [[ $seqs_to_keep == "keep_only_linked" ]]; then
    required_optional=$"required"
else
    required_optional=$"optional"
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
                echo "$fwd_primer;$required_optional...$rev_primer_RC;$required_optional" >> tempdir2/liked_fwd_revRC.fasta
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
                echo "$rev_primer;$required_optional...$fwd_primer_RC;$required_optional" >> tempdir2/liked_rev_fwdRC.fasta
                ((i=i+1))
            fi
        done < tempdir2/fwd_primer_RC.fasta
    fi
done < tempdir2/rev_primer.fasta

# make dirs for fwd and rev oriented seq batches
mkdir -p $output_dir/fwd_orient
mkdir -p $output_dir/rev_orient

### Read through each file in paired_end_files.txt
while read LINE; do
    #Write file name without extension
    inputR1=$(echo $LINE | sed -e "s/.$fileFormat//")
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
        untrimmed_output=$"--untrimmed-output $output_dir/untrimmed/$inputR1.$extension"
        untrimmed_paired_output=$"--untrimmed-paired-output $output_dir/untrimmed/$inputR2.$extension"
    fi

    ### Clip primers with cutadapt
    # when R1 starts with FWD primer
    if [[ $seqs_to_keep == "keep_all" ]]; then
        mkdir -p $output_dir/fwd_untrimmed
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        --untrimmed-output $output_dir/fwd_untrimmed/$inputR1.$extension \
        --untrimmed-paired-output $output_dir/fwd_untrimmed/$inputR2.$extension \
        --pair-filter=$pair_filter \
        -g file:tempdir2/fwd_primer.fasta \
        -G file:tempdir2/rev_primer.fasta \
        -o $output_dir/fwd_orient/$inputR1.$extension \
        -p $output_dir/fwd_orient/$inputR2.$extension \
        $inputR1.$extension $inputR2.$extension 2>&1)
        check_app_error

        # when R1 starts with REV primer
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=$pair_filter \
        -G file:tempdir2/fwd_primer.fasta \
        -g file:tempdir2/rev_primer.fasta \
        -o $output_dir/rev_orient/$inputR1.$extension \
        -p $output_dir/rev_orient/$inputR2.$extension \
        $output_dir/fwd_untrimmed/$inputR1.$extension $output_dir/fwd_untrimmed/$inputR2.$extension 2>&1)
        check_app_error

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        # when R1 starts with FWD primer
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        --untrimmed-output $output_dir/fwd_untrimmed/$inputR1.$extension \
        --untrimmed-paired-output $output_dir/fwd_untrimmed/$inputR2.$extension \
        --pair-filter=$pair_filter \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -G file:tempdir2/liked_rev_fwdRC.fasta \
        -o $output_dir/fwd_orient/$inputR1.$extension \
        -p $output_dir/fwd_orient/$inputR2.$extension \
        $inputR1.$extension $inputR2.$extension 2>&1)
        check_app_error

        # when R1 starts with REV primer
        checkerror=$(cutadapt --quiet \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=$pair_filter \
        -G file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/liked_rev_fwdRC.fasta \
        -o $output_dir/rev_orient/$inputR1.$extension \
        -p $output_dir/rev_orient/$inputR2.$extension \
        $output_dir/fwd_untrimmed/$inputR1.$extension $output_dir/fwd_untrimmed/$inputR2.$extension 2>&1)
        check_app_error
    fi
done < tempdir2/paired_end_files.txt

# delete empty files (if there are any)
find $output_dir/fwd_orient -empty -type f -delete
find $output_dir/rev_orient -empty -type f -delete

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"

if [[ -d "$output_dir/fwd_untrimmed" ]]; then
    rm -rf $output_dir/fwd_untrimmed
fi

### #stats for fwd_orient
cd $output_dir/fwd_orient 
#count reads before and after the process
touch seq_count_after.txt
outfile_check=$(ls *.$extension 2>/dev/null | wc -l)
if [[ $outfile_check != 0 ]]; then
    seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt
else
    printf '%s\n' "ERROR]: no fwd_orient output files generated ($output_dir). Not mixed orientation input?" >&2
    end_process
fi
touch /input/tempdir2/seq_count_before.txt
seqkit stats --threads 6 -T ../../*.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> /input/tempdir2/seq_count_before.txt
sed -i "s/\..\/\..\///" /input/tempdir2/seq_count_before.txt

#compile a track reads summary file (seq_count_summary.txt)
printf "File\tReads_in\tReads_out\n" > seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    while read LINE2; do
        file2=$(echo $LINE2 | awk '{print $1}')
        count2=$(echo $LINE2 | awk '{print $2}')
        if [[ "$file1" == "$file2" ]]; then
            printf "$file1\t$count1\t$count2\n" >> seq_count_summary.txt
        fi
    done < seq_count_after.txt
    #Report file where no sequences were reoriented (i.e. the output was 0)
    grep -Fq $file1 seq_count_after.txt
    if [[ $? != 0 ]]; then
        printf "$file1\t$count1\t0\n" >> seq_count_summary.txt
    fi
done < /input/tempdir2/seq_count_before.txt
rm seq_count_after.txt

### stats for rev_orient
cd $output_dir/rev_orient 
#count reads before and after the process
touch seq_count_after.txt
outfile_check=$(ls *.$extension 2>/dev/null | wc -l)
if [[ $outfile_check != 0 ]]; then
    seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt
else
    printf '%s\n' "ERROR]: no rev_orient output files generated ($output_dir). Not mixed orientation input?" >&2
    end_process
fi

#compile a track reads summary file (seq_count_summary.txt)
printf "File\tReads_in\tReads_out\n" > seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    while read LINE2; do
        file2=$(echo $LINE2 | awk '{print $1}')
        count2=$(echo $LINE2 | awk '{print $2}')
        if [[ "$file1" == "$file2" ]]; then
            printf "$file1\t$count1\t$count2\n" >> seq_count_summary.txt
        fi
    done < seq_count_after.txt
    #Report file where no sequences were reoriented (i.e. the output was 0)
    grep -Fq $file1 seq_count_after.txt
    if [[ $? != 0 ]]; then
        printf "$file1\t$count1\t0\n" >> seq_count_summary.txt
    fi
done < /input/tempdir2/seq_count_before.txt
rm seq_count_after.txt

### stats for untrimmed
cd $output_dir/untrimmed 
#count reads before and after the process
touch seq_count_after.txt
outfile_check=$(ls *.$extension 2>/dev/null | wc -l)
if [[ $outfile_check != 0 ]]; then
    seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt

    #compile a track reads summary file (seq_count_summary.txt)
    printf "File\tInput_reads\tUntrimmed_reads\n" > seq_count_summary.txt
    while read LINE; do
        file1=$(echo $LINE | awk '{print $1}')
        count1=$(echo $LINE | awk '{print $2}')
        while read LINE2; do
            file2=$(echo $LINE2 | awk '{print $1}')
            count2=$(echo $LINE2 | awk '{print $2}')
            if [[ "$file1" == "$file2" ]]; then
                printf "$file1\t$count1\t$count2\n" >> seq_count_summary.txt
            fi
        done < seq_count_after.txt
        #Report file where no sequences were reoriented (i.e. the output was 0)
        grep -Fq $file1 seq_count_after.txt
        if [[ $? != 0 ]]; then
            printf "$file1\t$count1\t0\n" >> seq_count_summary.txt
        fi
    done < /input/tempdir2/seq_count_before.txt
    rm seq_count_after.txt
fi

#Delete tempdir
if [[ -d "/input/tempdir2" ]]; then
    rm -rf /input/tempdir2
fi

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
printf "# Primers from mixed orient sequences were removed using cutadapt (see 'Core command' below for the used settings).

Files in 'primersCut_out/fwd_orient' folder represent forward orient sequences from where the PCR primers were recognized and clipped (R1 files started with forward primer and R2 files with reverse_complement of reverse primer).
Files in 'primersCut_out/rev_orient' folder represent reverse_complement sequences from where the PCR primers were recognized and clipped (R1 files started with reverse primer and R2 files with reverse_complement of forward primer).

Primers:
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].

Output R1 and R2 reads are synchronized for merging paired-end data. 

If no outputs were generated into $output_dir, check your specified primer strings and adjust settings.

Core commands -> \n" > $output_dir/README.txt

if [[ $seqs_to_keep == "keep_all" ]]; then
    printf "# when R1 starts with REV primer. seqs_to_keep == "keep_all"
cutadapt $mismatches $min_length $overlap $indels --untrimmed-output fwd_untrimmed/inputR1 --untrimmed-paired-output fwd_untrimmed/inputR2 --pair-filter=$pair_filter \
-g ForwardPrimer \
-G /ReversePrimer \
-o fwd_orient/inputR1. \
-p fwd_orient/inputR2 \
inputR1 inputR2 

# when R1 starts with REV primer
cutadapt $mismatches $min_length $overlap $indels $untrimmed_output $untrimmed_paired_output --pair-filter=$pair_filter \
-G ForwardPrimer \
-g ReversePrimer \
-o rev_orient/inputR1 \
-p rev_orient/inputR2 \
fwd_untrimmed/inputR1 fwd_untrimmed/inputR2 \n" >> $output_dir/README.txt

elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
    printf "# when R1 starts with FWD primer. seqs_to_keep == "keep_only_linked"
cutadapt $mismatches $min_length $overlap $indels --untrimmed-output $output_dir/fwd_untrimmed/$inputR1.$extension --untrimmed-paired-output $output_dir/fwd_untrimmed/$inputR2.$extension --pair-filter=$pair_filter \
-g liked_forwardPrimer_and_reverseComplementReversePrimer \
-G liked_reversePrimer_and_reverseComplementForwardPrimer \
-o fwd_orient/inputR1 \
-p fwd_orient/inputR2 \
inputR1 inputR2

# when R1 starts with REV primer
cutadapt $mismatches $min_length $overlap $indels $cores $untrimmed_output $untrimmed_paired_output --pair-filter=$pair_filter \
-G liked_forwardPrimer_and_reverseComplementReversePrimer \
-g liked_reversePrimer_and_reverseComplementForwardPrimer \
-o rev_orient/inputR1 \
-p rev_orient/inputR2 \
fwd_untrimmed/inputR1 fwd_untrimmed/inputR2. \n" >> $output_dir/README.txt
fi

printf "\nSummary of sequence counts in 'seq_count_summary.txt'

Total run time was $runtime sec.

##########################################################
###Third-party applications used for this process [PLEASE CITE]:
#cutadapt v4.4 for cutting the primers
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0 for generating reverse complementary primer strings
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" >> $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=paired_end"
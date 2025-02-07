#!/bin/bash

# REMOVE PRIMERS from paired-end reads, amplicons are expected to be MIXED oriented (5'-3' and 3'-5').  
# Degenerate primers are allowed using IUPAC codes.
# Input = paired-end fastq or paired-end fasta files. If using fasta, then cores must = 1
    # Using this only in the DADA2 paired-end pipeline where amplicons are MIXED oriented  - outputs the required structure for quality_filtering_paired_end_dada2_mixed.sh.
    # Output = files in primersCut_out/fwd_orient and primersCut_out/rev_orient dirs.

################################################
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
################################################
#Source for functions
source /scripts/submodules/framework.functions.sh

# Checking tool versions
cutadapt_version=$(cutadapt --version 2>&1)
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# Checking tool versions ...\n"
printf "# cutadapt (version $cutadapt_version)\n"
printf "# seqkit (version $seqkit_version)\n"

# load variables
extension=${fileFormat}  # KEEP THIS (removed in some other scripts)
mismatches=$"-e ${mismatches}"
min_length=$"--minimum-length 32"   # minimum len of the output sequence. FIXED to 32 (in order to avoid 0 len seqs) because cutadapt --minimum-length does not behave as expected
overlap=$"--overlap ${min_overlap}" # this option is ignored for anchored adapters since these do not allow partial matches.
cores=$"--cores ${cores}"
no_indels=$no_indels
discard_untrimmed=$"TRUE" #currently only fixed to TRUE 
seqs_to_keep=$seqs_to_keep
pair_filter=$pair_filter
fwd_tempprimer=$forward_primers
rev_tempprimer=$reverse_primers


#############################
### Start of the workflow ###
#############################
printf "# Running cut primers (from MIXED amplicons) \n"
if [[ $no_indels == "TRUE" ]]; then
    indels=$"--no-indels"
fi
# if keep_only_linked, then linked primers are REQUIRED (default = optional)
if [[ $seqs_to_keep == "keep_only_linked" ]]; then
    required_optional=$"required"
else
    required_optional=$"optional"
fi


# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
  echo "Working with multiple sequencing runs in multiRunDir"
  echo "Process = cut primers (from MIXED amplicons)"
  cd /input/multiRunDir
  # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")
    echo "Working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = cut primers (from MIXED amplicons)"
    DIRS=$(pwd)
fi

### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)
    cd $seqrun

    if [[ $multiDir == "TRUE" ]]; then
        ### Check if the dir has the specified file extension; if not then ERROR
        count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
        if [[ $count != 0 ]]; then
            printf "\n$ ---- Processing samples in $seqrun\n"
        else
            printf '%s\n' "ERROR]: cannot find files with specified extension '$fileFormat' in dir $seqrun.
            Please check the extension of your files and specify again.
            >Quitting" >&2
            end_process
        fi
        #output dir
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/primersCut_out"
        tempdir2=$"/input/multiRunDir/${seqrun%%/*}/tempdir2"
        input_files_path=$"/input/multiRunDir/${seqrun%%/*}"
        ### Prepare working env and check paired-end data
        prepare_PE_env
        # prepare primers
        prepare_primers
    else
        #output dir
        output_dir=$"/input/primersCut_out"
        tempdir2=$"/input/tempdir2"
        input_files_path=$"/input"
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check paired-end data
        prepare_PE_env
        # prepare primers
        prepare_primers
    fi

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
            -a file:tempdir2/liked_fwd_revRC.fasta \
            -G file:tempdir2/rev_primer.fasta \
            -A file:tempdir2/liked_rev_fwdRC.fasta \
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
            -g file:tempdir2/rev_primer.fasta \
            -a file:tempdir2/liked_rev_fwdRC.fasta \
            -G file:tempdir2/fwd_primer.fasta \
            -A file:tempdir2/liked_fwd_revRC.fasta \
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
    if (( $outfile_check != 0 )); then
        seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt
    else
        printf '%s\n' "ERROR]: no fwd_orient output files generated ($output_dir). Not mixed orientation input?" >&2
        end_process
    fi
    touch $tempdir2/seq_count_before.txt
    seqkit stats --threads 6 -T $input_files_path/*.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> $tempdir2/seq_count_before.txt

    #compile a track reads summary file (seq_count_summary.txt)
    printf "File\tReads_in\tReads_out\n" > seq_count_summary.txt
    while read LINE; do
        file1=$(basename $LINE) # get only the file name
        count1=$(echo $LINE | awk '{print $2}') # get the seq count
        while read LINE2; do
            file2=$(basename $LINE2)
            count2=$(echo $LINE2 | awk '{print $2}')
            if [[ "$file1" == "$file2" ]]; then
                printf "$file1\t$count1\t$count2\n" >> seq_count_summary.txt
            fi
        done < seq_count_after.txt
        #Report file where no sequences were outputted (i.e. the output was 0)
        grep -Fq $file1 seq_count_after.txt
        if [[ $? != 0 ]]; then
            printf "$file1\t$count1\t0\n" >> seq_count_summary.txt
        fi
    done < $tempdir2/seq_count_before.txt
    if [[ $debugger != "true" ]]; then
        rm seq_count_after.txt
    fi

    ### stats for rev_orient
    cd $output_dir/rev_orient 
    #count reads before and after the process
    touch seq_count_after.txt
    outfile_check=$(ls *.$extension 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then
        seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt
    else
        printf '%s\n' "ERROR]: no rev_orient output files generated ($output_dir). Not mixed orientation input?" >&2
        end_process
    fi

    #compile a track reads summary file (seq_count_summary.txt)
    printf "File\tReads_in\tReads_out\n" > seq_count_summary.txt
    while read LINE; do
        file1=$(basename $LINE)
        count1=$(echo $LINE | awk '{print $2}')
        while read LINE2; do
            file2=$(basename $LINE2)
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
    done < $tempdir2/seq_count_before.txt
    if [[ $debugger != "true" ]]; then
        rm seq_count_after.txt
    fi

    ### stats for untrimmed
    cd $output_dir/untrimmed 
    #count reads before and after the process
    touch seq_count_after.txt
    outfile_check=$(ls *.$extension 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then
        seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> seq_count_after.txt

        #compile a track reads summary file (seq_count_summary.txt)
        printf "File\tInput_reads\tUntrimmed_reads\n" > seq_count_summary.txt
        while read LINE; do
            file1=$(basename $LINE)
            count1=$(echo $LINE | awk '{print $2}')
            while read LINE2; do
                file2=$(basename $LINE2)
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
        done < $tempdir2/seq_count_before.txt
        if [[ $debugger != "true" ]]; then
            rm seq_count_after.txt
        fi
    fi
    # rm tempdir2
    if [[ $debugger != "true" ]]; then
        if [[ -d "$tempdir2" ]]; then
            rm -rf $tempdir2
        fi
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

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

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

################################################
###Third-party applications used for this process:
# cutadapt (version $cutadapt_version) for cutting the primers
   #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
   #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit (version $seqkit_version) for generating reverse complementary primer strings
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##############################################" >> $output_dir/README.txt

    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done 

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
else
    echo "workingDir=$output_dir"
fi
echo "fileFormat=$extension"
echo "readType=paired_end"

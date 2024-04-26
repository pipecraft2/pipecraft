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
extension=${fileFormat}             # KEEP THIS (removed in some other scripts)
mismatches=$"-e ${mismatches}"      # number of mismatches in primer search
min_length=$"--minimum-length 32"   # minimum len of the output sequence. FIXED to 32 (in order to avoid 0 len seqs) because cutadapt --minimum-length does not behave as expected
overlap=$"--overlap ${min_overlap}"
cores=$"--cores ${cores}"
# $no_indels                        # "--no-indels" disallows insertions and deletions in a primer
discard_untrimmed=$"TRUE"           # fixed to TRUE 
# $seqs_to_keep                     # all/keep_only_linked
# $pair_filter                      # both/any
fwd_tempprimer=$forward_primers
rev_tempprimer=$reverse_primers

#Source for functions
source /scripts/submodules/framework.functions.sh

#############################
### Start of the workflow ###
#############################
if [[ $no_indels == "true" ]]; then
    indels=$"--no-indels"
fi
# if keep_only_linked, then linked primers are REQUIRED (default = optional)
if [[ $seqs_to_keep == "keep_only_linked" ]]; then
    required_optional=$"required"
else
    required_optional=$"optional"
fi


# Check if I need to work with multiple or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
  echo "Working with multiple sequencing runs in multiRunDir"
  cd /input/multiRunDir
  # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir2" | grep -v "tempdir" | grep -v "skip_")
    echo "runs = $DIRS"
    multiDir=$"TRUE"
    export multiDir
else
  echo "Working with individual sequencing run"
  DIRS=$(pwd)
fi

### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start=$(date +%s)
    cd $seqrun
    
    if [[ $multiDir == "TRUE" ]]; then
        ### Check if the dir has the specified file extension; if not then reset the loop and check the next dir
        count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
        if [[ $count != 0 ]]; then
            printf "\n$ ---- Processing samples in $sequencing_run\n"
        else
            cd ..
            break 2
        fi
        #output dir
        output_dir=$"/input/multiRunDir/$seqrun/primersCut_out"
        ### Prepare working env and check paired-end data
        prepare_PE_env
        # prepare primers
        prepare_primers
    else
        #output dir
        output_dir=$"/input/primersCut_out"
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check paired-end data
        prepare_PE_env
        # prepare primers
        prepare_primers
    fi

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
            if [[ $debugger != "true" ]]; then
                rm $output_dir/$inputR1.round1.$extension $output_dir/$inputR2.round1.$extension
            fi

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
            if [[ $debugger != "true" ]]; then
                rm $output_dir/$inputR1.round1.$extension $output_dir/$inputR2.round1.$extension
            fi
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
    printf "# Primers were removed using cutadapt (see 'Core command' below for the used settings).

    Files in 'primersCut_out' folder represent sequences from where the PCR primers were recognized and clipped.
    Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
    Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
    [If primers were not specified in orientations noted above, please run this step again].

    Output R1 and R2 reads are synchronized for merging paired-end data. 
    If no outputs were generated into /$output_dir, check your specified primer strings and adjust settings.

    Core commnads -> \n" > $output_dir/README.txt

    if [[ $seqs_to_keep == "keep_all" ]]; then
        printf "seqs_to_keep == "keep_all"; cutadapt $mismatches $min_length $overlap $indels $untrimmed_output $untrimmed_paired_output --pair-filter=$pair_filter \
    -g liked_forwardPrimer_and_reverseComplementReversePrimer \
    -g liked_reversePrimer_and_reverseComplementForwardPrimer \
    -g ForwardPrimer \
    -a reverseComplementReversePrimer \
    -g ReversePrimer \
    -a reverseComplementForwardPrimer \
    -G liked_forwardPrimer_and_reverseComplementReversePrimer \
    -G liked_reversePrimer_and_reverseComplementForwardPrimer \
    -G ReversePrimer \
    -A reverseComplementForwardPrimer \
    -G forwardPrimer \
    -A reverseComplementReversePrimer \
    -o output.round1 -p output.round1 inputR1 inputR2 

    #round2; clipping if present, but no discarding
    cutadapt $mismatches $min_length $overlap $indels --pair-filter=$pair_filter \
    -g liked_forwardPrimer_and_reverseComplementReversePrimer \
    -g liked_reversePrimer_and_reverseComplementForwardPrimer \
    -g forwardPrimer \
    -a reverseComplementReversePrimer \
    -g ReversePrimer \
    -a reverseComplementForwardPrimer \
    -G liked_forwardPrimer_and_reverseComplementReversePrimer \
    -G liked_reversePrimer_and_reverseComplementForwardPrimer \
    -G ReversePrimer \
    -A reverseComplementForwardPrimer \
    -G forwardPrimer \
    -A reverseComplementReversePrimer \
    -o outputR1 -p outputR2 inputR1.from_round1 inputR2.from_round1 \n" >> $output_dir/README.txt

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        printf "seqs_to_keep == "keep_only_linked"; cutadapt $mismatches $min_length $overlap $indels $untrimmed_output $untrimmed_paired_output --pair-filter=$pair_filter \
    -g liked_forwardPrimer_and_reverseComplementReversePrimer \
    -g liked_reversePrimer_and_reverseComplementForwardPrimer \
    -G liked_forwardPrimer_and_reverseComplementReversePrimer \
    -G liked_reversePrimer_and_reverseComplementForwardPrimer \
    -o output.round1 -p output.round1 inputR1 inputR2 

    #additional check of the primer presence (no discarding)
    cutadapt $mismatches $min_length $overlap $indels \
    -g ForwardPrimer \
    -g reverseComplementForwardPrimer \
    -g ReversePrimer \
    -g reverseComplementReversePrimer \
    -G ForwardPrimer \
    -G reverseComplementForwardPrimer \
    -G ReversePrimer \
    -G reverseComplementReversePrimer \
    -o outputR1 -p outputR2 inputR1.from_round1 inputR2.from_round1 \n" >> $output_dir/README.txt
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

    
    ### if working with multiRunDir then cd ..
    if [[ $multiDir == "TRUE" ]]; then 
        cd ..
    fi
done 

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=paired_end"
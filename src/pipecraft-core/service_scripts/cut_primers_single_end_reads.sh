#!/bin/bash

# REMOVE PRIMERS from single-end reads
# Degenerate primers are allowed using IUPAC codes. Reverse complementary strings will be also searched.
# Input = single-end fastq or fasta files. If using fasta, then cores must = 1

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

if [[ $no_indels == "true" ]]; then
    indels=$"--no-indels"
fi
# if keep_only_linked, then linked primers are REQUIRED (default = optional)
if [[ $seqs_to_keep == "keep_only_linked" ]]; then
    required_optional=$"required"
else
    required_optional=$"optional"
fi

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/primersCut_out"

start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env

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
# Reverse complement REV primers
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

#############################
### Start of the workflow ###
#############################
### Read through each fastq/fasta file in folder
for file in *.$fileFormat; do
    #Write file name without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
    ## Preparing files
    printf "\n____________________________________\n"
    printf "Checking $input ...\n"

    ### Check input formats (only fastq/fasta supported)
    check_extension_fastxgz

    ##############################
    ### Start clipping primers ###
    ##############################
    printf "\n# Clipping primers from $input \n"
    printf " forward primer(s): $fwd_tempprimer\n"
    printf " reverse primer(s): $rev_tempprimer\n"

    #If discard_untrimmed = TRUE, then assigns outputs and make outdir
    if [[ $discard_untrimmed == "TRUE" ]]; then
        mkdir -p $output_dir/untrimmed
        untrimmed_output=$"--untrimmed-output $output_dir/untrimmed/$input.untrimmed.$extension"
    fi

    ### Clip primers with cutadapt
        # --revcomp compared RC seq strands, so no need to specify RC primers
    if [[ $seqs_to_keep == "keep_all" ]]; then
        checkerror=$(cutadapt --quiet --revcomp \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/fwd_primer.fasta \
        -a file:tempdir2/rev_primer_RC.fasta \
        -o $output_dir/$input.$extension \
        $input.$extension 2>&1)
        check_app_error

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        checkerror=$(cutadapt --quiet --revcomp \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -o $output_dir/$input.$extension \
        $input.$extension 2>&1)
        check_app_error
    fi
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
#file identifier string after the process
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for untrimmed seqs
if [[ $discard_untrimmed == "TRUE" ]]; then
    printf "Files in /untrimmed folder represent seqs that did not contain specified primer strings.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
If no files in this folder, then all sequences were passed to files in $output_dir directory" > $output_dir/untrimmed/README.txt
fi

#Make README.txt file for PrimerClipped reads
printf "# Primers were removed using cutadapt (see 'Core command' below for the used settings).

Files in $output_dir folder represent sequences from where the PCR primers were recognized and clipped.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].

Note that REVERSE COMPLEMENTARY search was also performed for sequences when no primer match was found on the 'original' strand ('--revcomp' setting).
If a match was found on a reverse complementary strand, then this reverse complementary sequence is outputted instead of 'original' read where no primer matches were found.
If forward primer(s) were specified in 5'-3' orientation, then all output seqs are in 5'-3' orientation.
Therefore, for single-end data, NO ADDITIONAL 'reorient reads' process is needed (and also impossible, because primers are now clipped).

If no outputs were generated into $output_dir, check your specified primer strings and adjust settings.

Core command -> \n" > $output_dir/README.txt

if [[ $seqs_to_keep == "keep_all" ]]; then
        printf "seqs_to_keep == "keep_all"; cutadapt --revcomp $mismatches $min_length $overlap $indels $untrimmed_output -g liked_forwardPrimer_and_reverseComplementReversePrimer -g forwardPrimer -a reverseComplementReversePrimer -o output input \n" >> $output_dir/README.txt
elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
    printf "seqs_to_keep == "keep_only_linked"; cutadapt --revcomp $mismatches $min_length $overlap $indels $untrimmed_output -g liked_forwardPrimer_and_reverseComplementReversePrimer -o output input \n" >> $output_dir/README.txt
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
echo "readType=single_end"

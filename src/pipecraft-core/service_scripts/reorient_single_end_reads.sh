#!/bin/bash

#Input = single-end fastq or fasta files.

#Reorient SINGLE-END reads: 
#mismatches = allowed number of differences for primer.
#Degenerate primers are allowed using IUPAC codes.

##########################################################
###Third-party applications:
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#fqgrep v0.4.4
    #Copyright (c) 2011-2016, Indraniel Das
    #https://github.com/indraniel/fqgrep
#pigz v2.4
##########################################################

###############################
#These variables are for testing (DELETE when implementing to PipeCraft)
extension=$fileFormat
mismatches=$mismatches
fwd_tempprimer=$forward_primers
rev_tempprimer=$reverse_primers
###############################
###############################

# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/reoriented_out"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_SE_env

### Process samples
for file in *.$extension; do
    ### Make temporary directory for temp files (for each sample)
    if [ -d tempdir ]; then
        rm -rf tempdir
    fi 
    mkdir tempdir
    #Write file name without extension
    input=$(echo $file | sed -e "s/.$extension//")
    outfile=$(echo $input | sed -e "s/.$extension//")
    ### Preparing files for reorienting
    printf "\n___________________________________\n"
    printf "Processing $file ...\n"

    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        #$extension will be $newextension
    check_gz_zip_SE
    ### Check input formats (fastx supported)
    check_extension_fastx

    #########################
    ### Start reorienting ###
    #########################
    printf "\n# Reorienting reads based on PCR primers: \n"
    ### Reorient FASTQ reads based on fwd and rev primers (up to 13 fwd and rev primers allowed). 
    SE_reorient_FWD &
    SE_reorient_REV & 
    wait
    #if rev primer found, then make reverse complementary and merge with 5_3.fastq file
    if [ -s tempdir/3_5.fastx ]; then
        checkerror=$(seqkit seq --quiet -t dna -r -p tempdir/3_5.fastx >> tempdir/5_3.fastx 2>&1)
        check_app_error
    fi

    ### Remove multiprimer artefacts
    multiprimer_search_SE

    # Move final output to $output_dir
    mkdir -p $output_dir
    mv tempdir/$outfile.$newextension $output_dir/$outfile.$newextension

    ### Move multiprimer chimeras to '$output_dir/multiprimer_chimeras' dir
    mkdir -p $output_dir/multiprimer_chimeras
    if [ -s tempdir/$input.multiprimer.$newextension ]; then
        mv tempdir/$input.multiprimer.$newextension $output_dir/multiprimer_chimeras
    fi

    ### Check if reoriented output is empty; if yes, then report WARNING
    if [ -s $output_dir/$outfile.$newextension ]; then
        :
    else
        printf '%s\n' "WARNING]: primers not found in file $outfile (no output)"
        rm $output_dir/$outfile.$newextension
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

#Make README.txt file for multi-primer chimeras
printf "If there are some files in that directory here,
then the sequences in these files were considered as chimeric ones and removed from the reoriented data set.\n\n
Logic behind considering these seqs as chimeric ones:
PCR primer strings were specified in orientation they are used in a PCR; i.e.
forward primer in 5'-3' orientation and reverse primer in 3'-5' orientation.
[IF THAT WAS NOT THE CASE, THEN RUN THIS STEP AGAIN!]
Therefore if a forward primer string was found in a sequence,
but also a reverse primer sting was found in the same sequence,
the sequence consists of 5'-3' and 3'-5' oriented fragments.
It is highly likely that this sequence is a chimeric one, and should be therefore removed.
Usually only very few such 'multi-primer' chimeric sequences are found in the amplicon data sets.\n" > $output_dir/multiprimer_chimeras/README.txt

#Make README.txt file for this process
printf "Files here represent sequences that have been reoriented based on PCR primers.
Forward primer(s) [has to be 5'-3']: $fwd_tempprimer
Reverse primer(s) [has to be 3'-5']: $rev_tempprimer
[If primers were not specified in orientations noted above, please run this step again].\n
RUNNING THE PROCESS SEVERAL TIMES IN THE SAME DIRECTORY WILL OVERWRITE ALL OUTPUTS!
\nSummary of sequence counts in 'seq_count_summary.txt'\n
\n\nTotal run time was $runtime sec.\n\n
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#seqkit v2.0.0 for manipulating reads
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
#fqgrep v0.4.4 for finding primer strings
    #Copyright (c) 2011-2016, Indraniel Das
    #https://github.com/indraniel/fqgrep
##########################################################" > $output_dir/README.txt

printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in 'seq_count_summary.txt'\n"
printf "Check README.txt file in $output_dir directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=/$output_dir"
echo "fileFormat=$newextension"
echo "dataFormat=demultiplexed"
echo "readType=single_end"

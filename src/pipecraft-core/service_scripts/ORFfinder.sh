#!/bin/bash

#Using ORFfinder to identify open reading frames.

regex='[^/]*$'
rep_seqs=${rep_seqs_file}
rep_seqs_file_path=$(echo $rep_seqs | grep -oP "$regex")
rep_seqs_temp=$(basename $rep_seqs_file_path) #basename, needed for macOS
rep_seqs_file=$(printf "/extraFiles/$rep_seqs_temp")
min_len=${min_len}
genetic_code=${genetic_code}
start_codon=${start_codon}
ignore_nested=${ignore_nested}
strand=${strand}

printf "Settings: 
rep_seqs_file=$rep_seqs_file
min_len=$min_len
genetic_code=$genetic_code
start_codon=$start_codon
ignore_nested=$ignore_nested
strand=$strand
"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)

#output dir
output_dir=$"/input/"

### Using ORFfinder to identify open reading frames
printf "Using ORFfinder to identify open reading frames ... \n"

ORFfinder -in $rep_seqs_file \
    -ml $min_len \
    -g $genetic_code \
    -s $start_codon \
    -n $ignore_nested \
    -strand $strand \
    -outfmt 1 \
    -out ORFfinder.temp

### retain the longest ORF
  # sort seqs by length
seqkit sort -l ORFfinder.temp \
    | sed -e 's/:.*//' | \
    sed -e 's/lcl|//' \
    > sort.temp
  # remove duplicate sequences (seqs with multiple ORFs). First seq will be kept
seqkit rmdup -n sort.temp -w 0 > ORFs.fasta

#write putative pseudogenes to a file
grep "^>" ORFs.fasta | sed -e 's/>//' > headers.temp
seqkit grep -w 0 -v -f headers.temp $rep_seqs_file > NUMTs.fasta

input_seqs=$(grep -c "^>" $rep_seqs_file)
numts=$(grep -c "^>" NUMTs.fasta)
orfs=$(grep -c "^>" ORFs.fasta)

printf "
input = $input_seqs seqs;
discarded $numts potential NUMTs;
output = $orfs seqs\n\n"

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
rm *.temp
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
printf "ORFfinder run for identifying open reading frames.

# ORFs.fasta = representative sequences from the input file with ORFs (open reading frames); contains $orfs reads.
# NUMTs.fasta = putative pseudogenes; contains $numts reads.

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#ORFfinder v0.4.3 for finding ORFs
    #https://www.ncbi.nlm.nih.gov/orffinder/
#seqkit v2.3.0 for sorting results
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" > $output_dir/README_ORFfinder.txt

#Done
printf "\nDONE"
printf "Data in directory $output_dir\n"
printf "Check README.txt file in $output_dir for further information about the process.\n"
printf "Total time: $runtime sec.\n"

#variables for all services
echo "workingDir=/$output_dir"
echo "fileFormat=$extension"
echo "dataFormat=demultiplexed"
echo "readType=single_end"



#!/bin/bash

#Using ORFfinder to identify open reading frames.
# Input = fasta file.

##################################
###Third-party applications:
#ORFfinder
#seqkit
#pigz
##################################

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

#output dir
output_dir=$"/input/"
# Source for functions
source /scripts/submodules/framework.functions.sh

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
# Check if rep_seqs is fasta
extension=$(echo $rep_seqs | (awk 'BEGIN{FS=OFS="."} {print $NF}';))
printf "\n fileFormat = $fileFormat \n"
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported here for 'rep seqs file'!
Supported extensions: fasta, fas, fa.
>Quitting" >&2
    end_process
fi

### Using ORFfinder to identify open reading frames
checkerror=$(ORFfinder -in $rep_seqs_file \
    -ml $min_len \
    -g $genetic_code \
    -s $start_codon \
    -n $ignore_nested \
    -strand $strand \
    -outfmt 1 \
    -out ORFfinder.temp 2>&1)
check_app_error

### retain the longest ORF
  # sort seqs by length
seqkit sort --quiet -l ORFfinder.temp \
    | sed -e 's/:.*//' | \
    sed -e 's/lcl|//' \
    > sort.temp
  # remove duplicate sequences (seqs with multiple ORFs). First seq will be kept
seqkit rmdup --quiet -n sort.temp -w 0 > ORFs.fasta

#write putative pseudogenes to a file
grep "^>" ORFs.fasta | sed -e 's/>//' > headers.temp
checkerror=$(seqkit grep --quiet -w 0 -v -f headers.temp $rep_seqs_file > NUMTs.fasta 2>&1)
check_app_error
checkerror=$(seqkit seq --quiet -n NUMTs.fasta > NUMTs_names.txt 2>&1)
check_app_error

# count outputs
input_seqs=$(grep -c "^>" $rep_seqs_file)
numts=$(grep -c "^>" NUMTs.fasta)
orfs=$(grep -c "^>" ORFs.fasta)

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
if [[ $debugger != "true" ]]; then
    rm *.temp
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
printf "Used ORFfinder to identify open reading frames (ORFs), kept the longest ORF sequence (see 'Core command' below for the used settings).

# ORFs.fasta = representative sequences from the input file with ORFs (open reading frames); contains $orfs reads.
# NUMTs.fasta = putative pseudogenes/off-target ORFs; contains $numts reads.
# NUMTs_names.txt = list of OTU/ASV identifiers of putative pseudogenes/off-target ORFs.

Core command -> 
ORFfinder -in $rep_seqs_file_path -ml $min_len -g $genetic_code -s $start_codon -n $ignore_nested -strand $strand
-ml = min_len
-g = genetic_code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
-s = start_codon
-n = ignore_nested

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#ORFfinder v0.4.3 for finding ORFs
    #https://www.ncbi.nlm.nih.gov/orffinder/
#seqkit v2.3.0 for sorting ORFfinder results
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" > $output_dir/README_ORFfinder.txt

#Done
printf "\nDONE"
printf "Total time: $runtime sec.\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
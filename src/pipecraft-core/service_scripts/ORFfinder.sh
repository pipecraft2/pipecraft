#!/bin/bash

#Using ORFfinder to identify open reading frames.
# Input = fasta file.

##################################
###Third-party applications:
#ORFfinder
#seqkit
#pigz
##################################
rep_seqs=${rep_seqs_file} # input fasta file
    regex='[^/]*$'
    rep_seqs_file_path=$(echo $rep_seqs | grep -oP "$regex")
    rep_seqs_temp=$(basename $rep_seqs_file_path) #basename, needed for macOS
    rep_seqs_file=$(printf "/extraFiles/$rep_seqs_temp")

# min_len       = positive integer. minimum length of the output sequence
# max_len       = positive integer. maximum length of the output sequence
# genetic_code  = positive integer. genetic code for translation. 5 = invertebrate mitochondrial code. Specify values from 1 to 33
# start_codon   = list [0, 1 or 2]. 0 = ATG only; 1 = ATG and alternative initation codons; 2 = any sense codon
# ignore_nested = bool. TRUE = ignore nested ORFs (completely placed within another)
# strand        = list [plus, minus, both]. output ORFs on specified strand only

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
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported here for 'rep seqs file'!
Supported extensions: fasta, fas, fa.
>Quitting" >&2
    end_process
fi

### Using ORFfinder to identify open reading frames
echo "# Running ORFfinder | "
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
echo "# Sorting results | "
in_name=$(echo $rep_seqs_file | awk -F "." '(NF=NF-1)' | sed -e 's/ /\./')
checkerror=$(seqkit sort --quiet -l ORFfinder.temp | sed -e 's/:.*//' | sed -e 's/lcl|//' > sort.temp 2>&1)
check_app_error
# remove duplicate sequences (seqs with multiple ORFs). First seq will be kept
checkerror=$(seqkit rmdup --quiet -n sort.temp | seqkit seq --remove-gaps --quiet -w 0 --max-len $max_len > $in_name.ORFs.fasta 2>&1)
check_app_error
checkerror=$(seqkit seq --quiet --name $in_name.ORFs.fasta > $in_name.ORFs.list 2>&1)
check_app_error


#write putative pseudogenes to a file
grep "^>" $in_name.ORFs.fasta | sed -e 's/>//' > headers.temp
checkerror=$(seqkit grep --quiet -w 0 -v -f headers.temp $rep_seqs_file > $in_name.NUMTs.fasta 2>&1)
check_app_error
checkerror=$(seqkit seq --quiet -n $in_name.NUMTs.fasta > $in_name.NUMTs_names.txt 2>&1)
check_app_error

# count outputs
input_seqs=$(grep -c "^>" $rep_seqs_file)
numts=$(grep -c "^>" $in_name.NUMTs.fasta)
orfs=$(grep -c "^>" $in_name.ORFs.fasta)

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
in=$(echo $in_name | sed -e 's/\/extraFiles\///')
printf "# Used ORFfinder to identify open reading frames (ORFs), kept the longest ORF per sequence (=<$max_len bp but =>$min_len bp) (see 'Core commands' below for the used settings).

Input file $in.$extension contained $input_seqs sequences.

Generated files:
# $in.ORFs.fasta      = representative sequences from the input file with ORFs (open reading frames); contains $orfs reads.
# $in.ORFs.list  = list of OTU/ASV ORF identifiers.
# $in.NUMTs.fasta     = putative pseudogenes/off-target ORFs; contains $numts reads.
# $in.NUMTs_names.txt = list of OTU/ASV identifiers of putative pseudogenes/off-target ORFs.

Core commands -> 
ORFfinder -in $rep_seqs_file_path -ml $min_len -g $genetic_code -s $start_codon -n $ignore_nested -strand $strand
seqkit seq -w 0 --max-len $max_len > $in.ORFs.fasta  # discard seqs >$max_len bp
    -ml = min_len
    -g = genetic_code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    -s = start_codon
    -n = ignore_nested
    --max-len = maximum length

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
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
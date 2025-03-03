#!/bin/bash

# Run seqkit stats on a folder with specified $fileFormat

##############################
### Third-party applications:
# vsearch
##############################
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# vsearch (version $vsearch_version)\n"

#get specified input fasta file
regex='[^/]*$'
fasta_file=$(echo $fasta_file | grep -oP "$regex")
fasta_file=$(printf "/extraFiles/$fasta_file")
printf "\n input fasta = $fasta_file \n"


# Source for functions
source /scripts/submodules/framework.functions.sh

#output dir
output_dir=$"/input/"

#############
### Start ###
#############
start_time=$(date)
start=$(date +%s)
fasta_basename=$(basename $fasta_file)

printf "# Running vsearch --allpairs_global \n"
checkerror=$(vsearch --allpairs_global $fasta_file \
        --id $id \
        --threads $cores \
        --userout $output_dir/$fasta_basename.pairwise_comparisons.txt \   
        --userfields query+target+id+alnlen+qcov+tcov 2>&1)
check_app_error



# 1. query: query label.
# 2. target: target (database sequence) label. The field is set to ’*’ if there is
# no alignment.
# 3. id: percentage of identity (real value ranging from 0.0 to 100.0). The percentage
# identity is defined as 100 * (matching columns) / (alignment
# length - terminal gaps). See fields id0 to id4 for other definitions.
# 4. alnlen: length of the query-target alignment (number of columns). The
# field is set to 0 if there is no alignment.
# 5. mism: number of mismatches in the alignment (zero or positive integer
# value).
# 6. opens: number of columns containing a gap opening (zero or positive integer
# value, excluding terminal gaps).
# 7. qlo: first nucleotide of the query aligned with the target. Always equal to
# 1 if there is an alignment, 0 otherwise (see qilo to ignore initial gaps).
# 8. qhi: last nucleotide of the query aligned with the target. Always equal to
# the length of the pairwise alignment, 0 otherwise (see qihi to ignore terminal
# gaps).
# 9. tlo: first nucleotide of the target aligned with the query. Always equal to
# 1 if there is an alignment, 0 otherwise (see tilo to ignore initial gaps).
# 10. thi: last nucleotide of the target aligned with the query. Always equal to
# the length of the pairwise alignment, 0 otherwise (see tihi to ignore terminal
# gaps).
# 11. evalue: expectancy-value (not computed for nucleotide alignments). Always
# set to -1.
# 12. bits: bit score (not computed for nucleotide alignments). Always set to 0.
# new userfields qilo, qihi, tilo, tihi give alignment coordinates ignoring terminal gaps,

### Add runtime and citation info to the output file
end=$(date +%s)
runtime=$((end-start))

printf "\n\n### This is tab-delimited file with the following columns:
    # file: name of the file
    # format: format of the file
    # type: type of the file (DNA/RNA)
    # num_seqs: number of sequences
    # sum_len: sum of the lengths of the sequences
    # min_len: minimum length of the sequences
    # avg_len: average length of the sequences
    # max_len: maximum length of the sequences

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

##############################################
###Third-party applications for this process:
# seqkit (version $seqkit_version)
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##############################################" >> $output_dir/seqkit_stats.$fileFormat.txt

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=$readType"
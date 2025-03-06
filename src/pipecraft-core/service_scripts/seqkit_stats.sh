#!/bin/bash

# Run seqkit stats on a folder with specified $fileFormat

##############################
### Third-party applications:
# seqkit
##############################
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# seqkit (version $seqkit_version)\n"
# Source for functions
source /scripts/submodules/framework.functions.sh

#output dir
output_dir=$"/input/"

#############
### Start ###
#############
start_time=$(date)
start=$(date +%s)
printf "# Running seqkit stats \n"
checkerror=$(seqkit stats --quiet -T --basename $output_dir/*.$fileFormat > $output_dir/seqkit_stats.$fileFormat.txt 2>&1)
check_app_error

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
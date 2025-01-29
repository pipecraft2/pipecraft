#!/bin/bash

# Merge sequencing runs processed with DADA2 if working with multuple runs in multiRunDir. 
 # Samples with the same name across runs are merged together.

##########################################################
###Third-party applications:
# dada2, R
##########################################################
start_time=$(date)
start=$(date +%s)
# source for functions
source /scripts/submodules/framework.functions.sh

# Excecute only if multiDir = true
if [[ ! -d "/input/multiRunDir" ]]; then
    printf '%s\n' "ERROR]: multiRunDir not detected. Cannot merge sequencing runs. Exiting.\n" >&2
    end_process
elif [[ $merge_runs == "true" ]]; then
    printf "Starting merge sequencing runs...\n"
    #output dir
    output_dir=$"/input/multiRunDir/merged_runs"
    export output_dir
    # remove output dir if it already exists
    if [[ -d "$output_dir" ]]; then
        rm -rf $output_dir
    fi
    # create new output dir
    mkdir -p $output_dir

    echo "input tables: $output_feature_table" 
    echo "input fasta: $output_fasta"
    echo "output dir: $output_dir"
else
    printf '%s\n' "ERROR]: Merge sequencing runs is not enabled. Exiting.\n" >&2
    end_process
fi

# Checking tool versions
R_version=$(R --version | head -n1 | cut -d " " -f3)
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# Checking tool versions ...\n"
printf "# R (version $R_version)\n"
printf "# DADA2 (version $dada2_version)\n"

# load variables
collapseNoMismatch=${collapseNoMismatch} # collapse identical ASVs (usearch_global --id 1)

### Merge ASV tables with dada2 mergeSequenceTables function in R
printf "# Running DADA2 mergeSequenceTables ...\n"
Rlog=$(Rscript /scripts/submodules/dada2_mergeRuns.R 2>&1)
echo $Rlog > $output_dir/dada2_mergeRuns.log 
wait
# format R log
sed -i 's/;; /\n/g' $output_dir/dada2_mergeRuns.log

# check if output files exist
if [[ ! -f $output_dir/ASVs.fasta ]] || [[ ! -f $output_dir/ASVs_table.txt ]]; then
    printf '%s\n' "ERROR]: Output files not found. Merge runs FAILED. 
    >Exiting." >&2
    end_process
fi

# count ASVs and sequences
ASV_count=$(grep -c "^>" $output_dir/ASVs.fasta)
nSeqs=$(awk 'BEGIN{FS=OFS="\t"}NR>1{for(i=2;i<=NF;i++) t+=$i; print t; t=0}' $output_dir/ASVs_table.txt  | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
nCols=$(awk -F'\t' '{print NF; exit}' $output_dir/ASVs_table.txt)
nSample=$(awk -v NUM=$nCols 'BEGIN {print (NUM-2)}') # -2 cuz 1st column is ASV_ID and 2nd is Sequence

### Make README.txt file (merged_runs)
end=$(date +%s)
runtime=$((end-start))
printf "# Merged sequencing runs with DADA2 mergeSequenceTables function.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input tables:\n" > $output_dir/README.txt

# Add each input table path
IFS=',' read -ra TABLES <<< "$output_feature_table"
for table in "${TABLES[@]}"; do
    printf "%s\n" "$table" >> $output_dir/README.txt
done

printf "
Output files:
------------
# ASVs_table.txt = merged ASV abundance table
# ASVs.fasta     = merged ASV sequences

Number of ASVs                       = $ASV_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample

###########################################################
###Third-party applications for this process [PLEASE CITE]:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
###########################################################" >> $output_dir/README.txt

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

# variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"

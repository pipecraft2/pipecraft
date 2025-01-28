#!/bin/bash

# Merge sequencing runs if working with multuple runs in multiRunDir. 
 # Samples with the same name across runs are merged together

##########################################################
###Third-party applications:
#vsearch, dada2, R
##########################################################

# source for functions
source /scripts/submodules/framework.functions.sh

# Excecute only if multiDir = true
if [[ ! -d "/input/multiRunDir" ]]; then
    printf '%s\n' "ERROR]: multiRunDir not detected. Cannot merge sequencing runs. Exiting.\n" >&2
    end_process
elif [[ $merge_runs == "true" ]]; then
    printf "Starting merge sequencing runs...\n" 
else
    printf '%s\n' "ERROR]: Merge sequencing runs is not enabled. Exiting.\n" >&2
    end_process
fi

echo "input tables: $output_feature_table" 
echo "input fasta: $output_fasta"

# Checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1)
seqkit_version=$(seqkit version 2>&1 | head -n 1)
R_version=$(R --version | head -n1 | cut -d " " -f3)
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | cut -d'"' -f2)
printf "# Checking tool versions ...\n"
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"
printf "# R (version $R_version)\n"
printf "# DADA2 (version $dada2_version)\n"

# load variables
collapseNoMismatch=${collapseNoMismatch} # collapse identical ASVs (usearch_global --id 1)

# collapse identical ASVs handling; for OTUs workflow [disabled in OTUs workflow]
if [[ $collapseNoMismatch == "" ]]; then
    collapseNoMismatch="false"
fi




#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"

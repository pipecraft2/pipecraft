#!/bin/bash

### VSEARCH SINTAX Taxonomy Assignment with Database Creation

################################################
###Third-party applications:
# vsearch
##############################################

# Checking tool versions
printf "# Checking tool versions ...\n"
vsearch_version=$(vsearch --version 2>&1 | head -n1 | awk '{print $2}')
printf "# VSEARCH version: %s\n" "$vsearch_version"

# Environment variables
workingDir=${workingDir}
extension=${fileFormat} && export fileFormat

# Load variables
fasta_file=${fasta_file}
sintax_db=${sintax_db}            # Provided SINTAX database file (FASTA or already built .udb)
sintax_cutoff=${sintax_cutoff}    # e.g., 0.8
sintax_strand=${sintax_strand}    # e.g., both
sintax_wordlength=${sintax_wordlength}  
sintax_threads=${sintax_threads}  

# Source for functions
source /scripts/submodules/framework.functions.sh

# Output directory
output_dir=$"/input/taxonomy_out.vsearch"
export output_dir

# Check and prepare SINTAX database
if [[ -f "$sintax_db" ]]; then
    ext="${sintax_db##*.}"
    if [[ "$ext" != "udb" ]]; then
        printf "# Building SINTAX database from FASTA file: %s\n" "$sintax_db"
        sintax_db_formatted="${sintax_db%.*}.udb"
        vsearch --makeudb_sintax "$sintax_db" --output "$sintax_db_formatted"
        if [[ $? -ne 0 ]]; then
            printf "Error: Failed to build SINTAX database from %s\n" "$sintax_db"
            exit 1
        fi
        sintax_db="$sintax_db_formatted"
    fi
else
    printf "Error: SINTAX database file %s does not exist.\n" "$sintax_db"
    exit 1
fi

# Start time
start_time=$(date)
start=$(date +%s)

### Check if files with specified extension exist in the directory
first_file_check
### Check if single-end files are compressed (decompress and check)
check_gz_zip_SE

#############################
### Start of the workflow ###
#############################
### Prepare working environment and check single-end data
prepare_SE_env

### Run VSEARCH SINTAX Taxonomy Assignment
printf "# Running VSEARCH SINTAX Taxonomy Assignment\n"
vsearch_log=$(vsearch --sintax "$fasta_file" \
        --db "$sintax_db" \
        --tabbedout "$output_dir/sintax_out.txt" \
        --sintax_cutoff "$sintax_cutoff" \
        --strand "$sintax_strand" \
        --wordlength "$sintax_wordlength" \
        --threads "$sintax_threads" 2>&1)
echo "$vsearch_log" > "$output_dir/vsearch_sintax.log"
printf "\nVSEARCH SINTAX Taxonomy Assignment completed\n"

########################################
### CLEAN UP AND COMPILE README FILE ###
########################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    if [[ -f $output_dir/vsearch_sintax.log ]]; then
        rm -f $output_dir/vsearch_sintax.log
    fi
fi

end=$(date +%s)
runtime=$((end - start))

### Make README.txt file
printf "# Taxonomy was assigned using VSEARCH SINTAX (with database preparation if required).
Query    = $fasta_file
Database = $sintax_db

# sintax_out.txt = VSEARCH SINTAX output file in tab-delimited format.
[If the provided database was not in the correct format, it was converted to a SINTAX database using 'vsearch --makeudb_sintax'.]

Core command -> 
vsearch --sintax $fasta_file --db $sintax_db --tabbedout sintax_out.txt --sintax_cutoff $sintax_cutoff --strand $sintax_strand --wordlength $sintax_wordlength --threads $sintax_threads

Total run time was $runtime sec.

################################################
###Third-party applications:
# vsearch (version $vsearch_version)
    #citation: Rognes, T., Flouri, T., Nichols, B. et al. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4: e2584. https://doi.org/10.7717/peerj.2584
##############################################" > $output_dir/README.txt

# Done
printf "\nDONE\n"
printf "Total time: $runtime sec.\n"

# Variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
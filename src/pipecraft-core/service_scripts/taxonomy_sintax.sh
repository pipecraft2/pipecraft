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
database_file=${database}      # Provided SINTAX database file (FASTA or already built .udb)
cutoff=${cutoff}          # e.g., default is 0.8
strand=${strand}          # [both, plus]
wordlength=${wordlength}  # Length of words (i.e. k-mers) for database indexing. Defaut is 8.
cores=${cores}            # number of cores to use

# input files (container paths)
database=$database_file
echo "database = $database"
echo "fasta_file = $fasta_file"
# overwrite fileFormat variable; get it from input fasta_file
fileFormat=$(echo $fasta_file | awk -F. '{print $NF}')
export fileFormat

# Source for functions
source /scripts/submodules/framework.functions.sh

# Output directory
output_dir=$"/input/taxonomy_out.sintax"
export output_dir

# Check and prepare SINTAX database
# The database file is mounted read-only; a .udb built from a fasta file is
# kept in PipeCraft's database cache (submodules/db_cache.sh).
source /scripts/submodules/db_cache.sh
if [[ ! -f "$database" ]]; then
    printf "Error: SINTAX database file %s does not exist.\n" "$database"
    exit 1
fi
ext=$(echo "${database##*.}" | tr '[:upper:]' '[:lower:]')
if [[ "$ext" != "udb" ]]; then
    db_cache=$(db_cache_dir)
    if db_cache_is_complete "$db_cache"; then
        printf "# Using cached SINTAX database for %s\n" "$(basename "$database")"
    else
        printf "# Building SINTAX database from FASTA file: %s (saved for future runs)\n" "$(basename "$database")"
        db_cache_reset "$db_cache"
        vsearch --makeudb_usearch "$database" --output "$db_cache/db.udb"
        if [[ $? -ne 0 ]]; then
            printf "Error: Failed to build SINTAX database from %s\n" "$database"
            exit 1
        fi
        db_cache_mark_complete "$db_cache"
    fi
    database_name=$(basename "$database")
    database="$db_cache/db.udb"
fi

# Start time
start_time=$(date)
start=$(date +%s)

#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
check_gz_zip_SE
### Check input formats (fasta supported)
check_extension_fasta

#############################
### Start of the workflow ###
#############################
echo "output_dir = $output_dir"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir

### Run VSEARCH SINTAX Taxonomy Assignment
printf "# Running VSEARCH SINTAX \n"
checkerror=$(vsearch --sintax "$fasta_file" \
        --db "$database" \
        --tabbedout "$output_dir/taxonomy.sintax.txt" \
        --sintax_cutoff "$cutoff" \
        --strand "$strand" \
        --wordlength "$wordlength" \
        --threads "$cores" 2>&1)
check_app_error
printf "\n SINTAX completed\n"

########################################
### CLEAN UP AND COMPILE README FILE ###
########################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
fi

end=$(date +%s)
runtime=$((end - start))

### Make README.txt file
# File names only, for the README
fasta_file=$(basename "$fasta_file")
database=${database_name:-$(basename "$database")}

printf "# Taxonomy was assigned using VSEARCH SINTAX (with database preparation if required).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Query    = $fasta_file
Database = $database

# taxonomy.sintax.txt = VSEARCH SINTAX output file in tab-delimited format.
[If the provided database was not in the correct format, it was converted to a SINTAX database using 'vsearch --makeudb_usearch'.]

Core command -> 
vsearch --sintax $fasta_file --db $database --tabbedout taxonomy.sintax.txt --sintax_cutoff $cutoff --strand $strand --wordlength $wordlength --threads $cores

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
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
sintax_cutoff=${sintax_cutoff}    
sintax_strand=${sintax_strand}    
sintax_wordlength=${sintax_wordlength}  
sintax_threads=${sintax_threads}  

# Source for functions
source /scripts/submodules/framework.functions.sh

# Function to check if database is in SINTAX format
check_sintax_format() {
    local db_path="$1"
    
    # If it's already a UDB file, it's valid
    if [[ "${db_path##*.}" == "udb" ]]; then
        return 0
    fi
    
    # Check if file exists
    if [[ ! -f "$db_path" ]]; then
        printf "Error: Database file '%s' does not exist.\n" "$db_path"
        return 1
    fi
    
    # Check if file is FASTA format 
    if [[ $(head -c 1 "$db_path") != ">" ]]; then
        printf "Error: Database does not appear to be in FASTA format (should start with '>').\n"
        return 1
    fi
    
    # Get only the first 10 headers
    local first_ten_headers=$(head -n 1000 "$db_path" | grep "^>" | head -n 10)
    
    # Check file size to decide if we need to check last headers
    local file_size=$(stat -c %s "$db_path")
    
    # Initialize counts
    local valid_count=0
    local invalid_count=0
    local first_invalid=""
    
    # Process first 10 headers
    while IFS= read -r header; do
        if [[ "$header" == *";tax="* ]] && [[ "$header" =~ \;tax=.*[dpcosfg]: ]]; then
            valid_count=$((valid_count + 1))
        else
            invalid_count=$((invalid_count + 1))
            if [[ -z "$first_invalid" ]]; then
                first_invalid="$header"
            fi
        fi
    done <<< "$first_ten_headers"
    
    # If file is large enough, check last headers too (only if needed)
    if [[ "$file_size" -gt 50000 && "$invalid_count" -eq 0 ]]; then
        # Read last ~10 headers
        local last_ten_headers=$(tail -n 2000 "$db_path" | grep "^>" | tail -n 10)
        
        # Check last headers only if all first headers were valid
        while IFS= read -r header; do
            if [[ "$header" == *";tax="* ]] && [[ "$header" =~ \;tax=.*[dpcosfg]: ]]; then
                valid_count=$((valid_count + 1))
            else
                invalid_count=$((invalid_count + 1))
                if [[ -z "$first_invalid" ]]; then
                    first_invalid="$header"
                fi
            fi
        done <<< "$last_ten_headers"
    fi
    
    # Determine if database is valid
    if [[ "$valid_count" -gt 0 && "$invalid_count" -eq 0 ]]; then
        return 0
    else
        printf "Error: Database is not in SINTAX format.\n"
        if [[ -n "$first_invalid" ]]; then
            printf "Example of invalid header: %s\n" "$first_invalid"
        fi
        printf "Headers should contain ';tax=d:Domain,p:Phylum,...;'\n"
        return 1
    fi
}

# Output directory
output_dir=$"/input/taxonomy_out.vsearch"
export output_dir

# Check and prepare SINTAX database
if [[ -f "$sintax_db" ]]; then
    ext="${sintax_db##*.}"
    if [[ "$ext" != "udb" ]]; then
        printf "# Building SINTAX database from FASTA file: %s\n" "$sintax_db"
        sintax_db_formatted="${sintax_db%.*}.udb"
        vsearch --makeudb_usearch "$sintax_db" --output "$sintax_db_formatted"
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
check_app_error

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
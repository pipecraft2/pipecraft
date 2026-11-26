#!/bin/bash
set -euo pipefail

# Inputs
#   table      = feature table file
#   fasta_file = corresponding FASTA file

: "${table:?ERROR: missing env var 'table'}"
: "${fasta_file:?ERROR: missing env var 'fasta_file'}"

# Source for functions
source /scripts/submodules/framework.functions.sh

# Inputs (container paths)
table_in=$table
export table_in
printf "\n input table = %s \n" "$table_in"

fasta_in=$fasta_file
export fasta_in
printf "\n input fasta = %s \n" "$fasta_in"

start_time=$(date)
start=$(date +%s)

# Run R script
Rscript /scripts/submodules/add_sequences_to_table.R --table "$table_in" --fasta_file "$fasta_in"

# Calculate runtime
end=$(date +%s)
runtime=$((end - start))

# Output directory (PipeCraft2 convention)
output_dir="/input"
export output_dir

# Derive output table name (matches the R script logic)
base_name=$(basename "$table")
base_name="${base_name%.*}"
output_name="${base_name}_wSeqs.txt"

# Make README.txt file (PipeCraft2 style)
printf "# Sequences were added to the feature table (2nd column) from the corresponding FASTA file.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input table  = $(basename "$table_in")
Input FASTA  = $(basename "$fasta_in")
Output table = $output_name
" > "$output_dir/README_add_sequences_to_table.txt"



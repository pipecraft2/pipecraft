#!/bin/bash

#Input = fasta file with sequences to query against BOLD Systems v5.
#Output = taxonomy assignment results in taxonomy_out.boldigger3/

##########################################################
###Third-party applications:
#boldigger3
    #citation: Buchner D, Leese F (2020) BOLDigger – a Python package to identify
    #  and organise sequences with the Barcode of Life Data systems.
    #  Metabarcoding and Metagenomics 4: e53535.
    #  https://doi.org/10.3897/mbmg.4.53535
##########################################################

# Check tool version
printf "# Checking tool versions ...\n"
boldigger3_version=$(pip show boldigger3 2>/dev/null | grep "^Version" | awk '{print $2}')
printf "# boldigger3 version: %s\n" "$boldigger3_version"

# Environment variables
workingDir=${workingDir}
extension=${fileFormat} && export fileFormat

# Load variables
fasta_file=${fasta_file}
database=${database}      # BOLD database number (1-8)
mode=${mode}              # Operating mode (1-3)
thresholds=${thresholds}  # Space-separated thresholds (optional)

# Prep input fasta file path for container
regex='[^/]*$'
fasta_basename=$(echo $fasta_file | grep -oP "$regex")
fasta_in=$(printf "/extraFiles/$fasta_basename")
echo "fasta_file = $fasta_in"
fasta_base=$(basename "$fasta_in" .fasta)

# Source for functions
source /scripts/submodules/framework.functions.sh

# Output directory
output_dir=$"/input/taxonomy_out.boldigger3"
export output_dir

#############################
### Start of the workflow ###
#############################
printf "output_dir = %s\n" "$output_dir"

# Remove old results and any stale intermediate files from previous runs
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
rm -rf /input/boldigger3_data 2>/dev/null

mkdir -p $output_dir

# Start time
start_time=$(date)
start=$(date +%s)

# Build thresholds flag (optional)
if [[ -n "$thresholds" ]] && [[ "$thresholds" != "undefined" ]] && [[ "$thresholds" != "null" ]]; then
    thresholds_flag="--thresholds $thresholds"
else
    thresholds_flag=""
fi

# Copy fasta into output_dir so BOLDigger3 creates boldigger3_data/ inside output_dir
# (BOLDigger3 always writes output relative to the parent directory of the input fasta path)
fasta_local="$output_dir/$fasta_base.fasta"
cp "$fasta_in" "$fasta_local"

### Run BOLDigger3
printf "# Running BOLDigger3 identify ...\n"
checkerror=$(boldigger3 identify "$fasta_local" \
    --db "$database" \
    --mode "$mode" \
    $thresholds_flag 2>&1)
check_app_error

printf "\n BOLDigger3 completed\n"

########################################
### COMPILE RESULTS AND CLEAN UP     ###
########################################

# Move all BOLDigger3 result files from boldigger3_data/ into output_dir
boldigger3_data_dir="$output_dir/boldigger3_data"

if [[ -d "$boldigger3_data_dir" ]]; then
    # Keep full output parity with standalone BOLDigger3 runs.
    cp -a "$boldigger3_data_dir"/. "$output_dir/"
    rm -rf "$boldigger3_data_dir"
fi

# Remove the temporary fasta copy
rm -f "$fasta_local"

end=$(date +%s)
runtime=$((end-start))

# Make README.txt file
printf "# Taxonomy was assigned using BOLDigger3 (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Query      = $fasta_basename
Database   = $database
Mode       = $mode
Thresholds = ${thresholds:-default (97 95 90 85)}

${fasta_base}_identification_result.xlsx           = BOLDigger3 top-hit taxonomy assignments (Excel format)
${fasta_base}_identification_result.parquet.snappy = BOLDigger3 top-hit taxonomy assignments (Parquet format)
Additional BOLDigger3-generated files (e.g. metadata/cache/intermediate files) are also retained.

Databases:
  1 = ANIMAL LIBRARY (PUBLIC)
  2 = ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)
  3 = ANIMAL LIBRARY (PUBLIC + PRIVATE)
  4 = VALIDATED CANADIAN ARTHROPOD LIBRARY
  5 = PLANT LIBRARY (PUBLIC)
  6 = FUNGI LIBRARY (PUBLIC)
  7 = ANIMAL SECONDARY MARKERS (PUBLIC)
  8 = VALIDATED ANIMAL RED LIST LIBRARY

Operating modes:
  1 = Rapid Species Search  (up to 1000 seqs/batch, ~10000 seqs/hour)
  2 = Genus and Species Search (200 seqs/batch)
  3 = Exhaustive Search (100 seqs/batch, most thorough)

Default similarity thresholds: Species=97%%, Genus=95%%, Family=90%%, Order=85%%
  Custom thresholds: pass up to 5 space-separated values for Species Genus Family Order Class.

Core command ->
boldigger3 identify $fasta_basename --db $database --mode $mode $thresholds_flag

##########################################################
###Third-party applications [PLEASE CITE]:
#boldigger3 (version $boldigger3_version)
    #citation: Buchner D, Leese F (2020) BOLDigger – a Python package to identify
    #  and organise sequences with the Barcode of Life Data systems.
    #  Metabarcoding and Metagenomics 4: e53535.
    #  https://doi.org/10.3897/mbmg.4.53535
##########################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

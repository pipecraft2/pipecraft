#!/bin/bash

### Taxonomy assignment with BOLDigger3
### Queries sequences against BOLD Systems v5 online databases.
### Outputs: .xlsx and .parquet result files in taxonomy_out.boldigger3/

##########################################################
###Third-party applications:
# boldigger3
    #citation: Buchner D, Leese F (2020) BOLDigger – a Python package to identify
    #  and organise sequences with the Barcode of Life Data systems.
    #  Metabarcoding and Metagenomics 4: e53535.
    #  https://doi.org/10.3897/mbmg.4.53535
##########################################################

# Clean up any stale containers from previous runs to prevent naming conflicts
docker container rm $(docker ps -a --filter "name=BOLDigger3" -q) 2>/dev/null || true

# Check tool version
printf "# Checking tool versions ...\n"
boldigger3_version=$(pip show boldigger3 2>/dev/null | grep "^Version" | awk '{print $2}')
printf "# BOLDigger3 version: %s\n" "$boldigger3_version"

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

# Source for functions
source /scripts/submodules/framework.functions.sh

# Output directory
output_dir=$"/input/taxonomy_out.boldigger3"
export output_dir

#############################
### Start of the workflow ###
#############################
printf "output_dir = %s\n" "$output_dir"

# Remove old results if they exist to ensure clean new run
if [[ -d $output_dir ]]; then
    printf "# Removing old results from previous run ...\n"
    rm -rf $output_dir
fi

# Also clean up any stale intermediate files from previous runs in /input
printf "# Cleaning up intermediate files from previous runs ...\n"
rm -f /input/*.duckdb 2>/dev/null || true
rm -f /input/*.pkl 2>/dev/null || true
rm -f /input/*_bold_results*.xlsx 2>/dev/null || true
rm -f /input/*_identification_result.* 2>/dev/null || true

# Create fresh output directory
mkdir -p $output_dir
printf "# Output directory ready: $output_dir\n"

# Start time
start_time=$(date)
start=$(date +%s)

# Build thresholds flag (optional)
if [[ -n "$thresholds" ]] && [[ "$thresholds" != "undefined" ]] && [[ "$thresholds" != "null" ]]; then
    thresholds_flag="--thresholds $thresholds"
else
    thresholds_flag=""
fi

### Run BOLDigger3 from the output directory so cache and results land there
printf "# Running BOLDigger3 identify ...\n"
cd $output_dir

checkerror=$(boldigger3 identify "$fasta_in" \
    --db "$database" \
    --mode "$mode" \
    $thresholds_flag 2>&1)
check_app_error

printf "\n BOLDigger3 completed\n"

########################################
### CONSOLIDATE AND CLEAN UP OUTPUT ###
########################################

# BOLDigger3 creates files with basename from input FASTA
# Extract base name without extension for pattern matching (e.g., "test" from "test.fasta")
fasta_base=$(basename "$fasta_in" .fasta)

# Find the actual result files that BOLDigger3 generated
# BOLDigger3 outputs: *_identification_result.xlsx and *_identification_result.parquet.snappy
find /input -maxdepth 1 -name "*_identification_result.*" 2>/dev/null | while read result_file; do
    if [[ -f "$result_file" ]]; then
        # Extract just the filename
        result_basename=$(basename "$result_file")
        # Copy to output directory, renaming to a clean name based on fasta basename
        filename_nobase=$(echo "$result_basename" | sed "s/^[^_]*_//")  # Remove everything before first underscore
        cp "$result_file" "$output_dir/$fasta_base.$filename_nobase"
        printf "  Moved: $result_basename → $fasta_base.$filename_nobase\n"
    fi
done

# Clean up intermediate/cache files that are not needed
printf "# Cleaning up intermediate files ...\n"
rm -f /input/*.duckdb 2>/dev/null                           # Remove DuckDB cache
rm -f /input/*_bold_results*.xlsx 2>/dev/null              # Remove intermediate BOLD xlsx files
rm -f /input/*.pkl 2>/dev/null                             # Remove pickle cache files
find "$output_dir" -name "*_bold_results*.xlsx" -delete 2>/dev/null  # Remove intermediate files in output
find "$output_dir" -name "*.duckdb" -delete 2>/dev/null    # Remove duckdb if any in output

end=$(date +%s)
runtime=$((end - start))

# Generate clean output file list
identification_xlsx=$(ls -1 "$output_dir"/*identification_result.xlsx 2>/dev/null | head -1)
identification_parquet=$(ls -1 "$output_dir"/*identification_result.parquet.snappy 2>/dev/null | head -1)

### Make comprehensive README.txt file
printf "# Taxonomy Assignment using BOLDigger3 (v${boldigger3_version})
# Query against BOLD Systems v5 online database
#
# Execution Info:
Start time:  $start_time
End time:    $(date)
Total time:  $runtime seconds

# Input Parameters:
Query file:    $fasta_basename
Database:      $database
Mode:          $mode
Thresholds:    ${thresholds:-default (97 95 90 85)}

## OUTPUT FILES

### Primary Result Files:
$(basename "$identification_xlsx" 2>/dev/null || echo "[identification_result.xlsx]")
  - Excel format containing BOLD top-hit taxonomy assignments
  - Columns: sequence IDs, taxonomy (Phylum/Class/Order/Family/Genus/Species)
  - Includes BOLDigger3 flagging column for uncertainty indicators

$(basename "$identification_parquet" 2>/dev/null || echo "[identification_result.parquet.snappy]")
  - Parquet format (snappy-compressed) of the same data
  - Use for downstream bioinformatics pipelines or data science workflows
  - More efficient than Excel for large datasets

## REFERENCE INFORMATION

### Databases (--db parameter):
  1 = ANIMAL LIBRARY (PUBLIC)
  2 = ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)
  3 = ANIMAL LIBRARY (PUBLIC + PRIVATE)
  4 = VALIDATED CANADIAN ARTHROPOD LIBRARY
  5 = PLANT LIBRARY (PUBLIC)
  6 = FUNGI LIBRARY (PUBLIC)
  7 = ANIMAL SECONDARY MARKERS (PUBLIC)
  8 = VALIDATED ANIMAL RED LIST LIBRARY

### Operating Modes (--mode parameter):
  1 = Rapid Species Search  (up to 1000 seqs/batch, ~10000 seqs/hour) - fastest
  2 = Genus and Species Search (200 seqs/batch) - moderate speed/thoroughness
  3 = Exhaustive Search (100 seqs/batch) - slowest, most thorough

### Default Similarity Thresholds:
  Species: 97%, Genus: 95%, Family: 90%, Order: 85%

### BOLDigger3 Flagging System:
  Flag 1: Reverse BIN Taxonomy (species assigned via reverse BIN)
  Flag 2: Differing Taxonomic Information (<90% hit agreement)
  Flag 3: Private Data (all top 100 hits are private)
  Flag 4: Unique Hit (single hit among top 100)
  Flag 5: Multiple BINs (species-level hit spans multiple BINs)

## COMMAND USED
boldigger3 identify $fasta_basename --db $database --mode $mode $thresholds_flag

## CITATION
Buchner D, Leese F (2020) BOLDigger – a Python package to identify
and organise sequences with the Barcode of Life Data systems.
Metabarcoding and Metagenomics 4: e53535.
https://doi.org/10.3897/mbmg.4.53535

GitHub: https://github.com/DominikBuchner/BOLDigger3
BOLD Systems: https://www.boldsystems.org/" > $output_dir/README.txt

# Done
printf "\nDONE\n"
printf "Total time: $runtime sec.\n"

# Variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

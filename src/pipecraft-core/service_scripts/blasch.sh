#!/bin/bash

##BlasCh for chimera detection and recovery
##input: working directory with *.chimeras.fasta files and sample FASTA files
##output: rescued sequences and analysis reports

#load variables
extension=$fileFormat
dataFormat=$dataFormat
readType=$readType
workingDir=$workingDir
fileFormat=$fileFormat
blasch_output_dir=$workingDir/BlasCh_out

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dirs
mkdir -p $blasch_output_dir

#############################
#  BLAST-based Chimera      #
#  Detection and Recovery   #
#############################

printf "\n\nBlasCh: False positive chimera detection and recovery\n"
printf "Working directory: $workingDir\n"
printf "Output directory: $blasch_output_dir\n"

# Use working directory for input chimeras
input_chimeras=$workingDir

# Check if .chimeras.fasta files exist in working directory
chimeras_count=$(find "$input_chimeras" -name "*.chimeras.fasta" -type f | wc -l)
if [[ $chimeras_count -eq 0 ]]; then
    printf "\nERROR: No .chimeras.fasta files found in working directory: $input_chimeras\n"
    printf "Please ensure that chimera detection has been run prior to BlasCh analysis.\n"
    end_process
fi

# Use working directory for self FASTA files (auto-create from .fasta files)
self_fasta=$workingDir

# Check for Python script availability
# First try in the local directory structure
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
BLASCH_SCRIPT="${SCRIPT_DIR}/submodules/Blasch_PipeCraft.py"

if [[ -f "${BLASCH_SCRIPT}" ]]; then
    printf "Found BlasCh Python script at ${BLASCH_SCRIPT}\n"
elif [[ -f "/scripts/submodules/Blasch_PipeCraft.py" ]]; then
    # Fall back to Docker container path
    BLASCH_SCRIPT="/scripts/submodules/Blasch_PipeCraft.py"
    printf "Found BlasCh Python script at ${BLASCH_SCRIPT}\n"
else
    printf "\nERROR: BlasCh Python script not found at either:\n"
    printf "  - ${SCRIPT_DIR}/submodules/Blasch_PipeCraft.py\n"
    printf "  - /scripts/submodules/Blasch_PipeCraft.py\n"
    end_process
fi

# Set up BlasCh command arguments
blasch_args="--input_chimeras_dir $input_chimeras"
blasch_args="$blasch_args --self_fasta_dir $self_fasta"
blasch_args="$blasch_args --output_dir $blasch_output_dir"
blasch_args="$blasch_args --threads $threads"

# Add threshold parameters
blasch_args="$blasch_args --high_coverage_threshold $coverage"
blasch_args="$blasch_args --high_identity_threshold $identity"
blasch_args="$blasch_args --borderline_coverage_threshold $borderline_coverage"
blasch_args="$blasch_args --borderline_identity_threshold $borderline_identity"

# Add reference database if provided
if [[ "$reference_db" != "undefined" ]] && [[ -n "$reference_db" ]]; then
    # Convert to container path format
    reference_db_name=$(basename "$reference_db")
    container_reference_db="/extraFiles/$reference_db_name"
    
    # Check if the file exists in the container mount point
    if [[ -f "$container_reference_db" ]]; then
        blasch_args="$blasch_args --reference_db $container_reference_db"
        printf "Using reference database: $container_reference_db\n"
    else
        printf "WARNING: reference_db file not found in container: $container_reference_db. Proceeding without reference database.\n"
        printf "Original path was: $reference_db\n"
    fi
else
    printf "No reference database specified. Using only self-databases.\n"
fi

printf "\nRunning BlasCh with the following parameters:\n"
printf "Input chimeras directory: $input_chimeras\n"
printf "Self FASTA directory: $self_fasta\n"
printf "Output directory: $blasch_output_dir\n"
printf "Threads: $threads\n"
printf "High coverage threshold: $coverage%%\n"
printf "High identity threshold: $identity%%\n"
printf "Borderline coverage threshold: $borderline_coverage%%\n"
printf "Borderline identity threshold: $borderline_identity%%\n"

# Run BlasCh
printf "\nExecuting Python script: $BLASCH_SCRIPT\n"
printf "Full command: python3 $BLASCH_SCRIPT $blasch_args\n"
checkerror=$(python3 "$BLASCH_SCRIPT" $blasch_args 2>&1)
check_app_error

printf "\nBlasCh analysis completed successfully!\n"
printf "Results are available in: $blasch_output_dir\n"
printf "\nOutput files:\n"
printf "- Non-chimeric sequences: non_chimeric/*_non_chimeric.fasta\n"
printf "- Borderline sequences: borderline/*_borderline.fasta\n"
printf "- Chimeric sequences: detailed_results/*_chimeric.fasta\n"
printf "- Multiple alignment sequences: detailed_results/*_multiple_alignments.fasta\n"
printf "- Analysis report: chimera_recovery_report.txt\n"
printf "- Detailed results: detailed_results/*_sequence_details.csv\n"
printf "- Compressed XML files: xml/blast_results.zip\n"

#end
printf "\nDone\n"

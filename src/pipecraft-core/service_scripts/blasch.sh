#!/bin/bash

##BlasCh for chimera detection and recovery
##input: directory with *.chimeras.fasta files
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

# Check if input chimeras directory is provided and exists
if [[ "$input_chimeras_dir" == "undefined" ]]; then
    printf "\nERROR: input_chimeras_dir not specified\n"
    end_process
fi

if [[ ! -d "$input_chimeras_dir" ]]; then
    printf "\nERROR: input_chimeras_dir does not exist: $input_chimeras_dir\n"
    end_process
fi

# Check if self_fasta_dir is provided and exists
if [[ "$self_fasta_dir" == "undefined" ]]; then
    printf "\nERROR: self_fasta_dir not specified\n"
    end_process
fi

if [[ ! -d "$self_fasta_dir" ]]; then
    printf "\nERROR: self_fasta_dir does not exist: $self_fasta_dir\n"
    end_process
fi

# Check for Python script availability
if [[ ! -f "/scripts/submodules/Blasch_PipeCraft.py" ]]; then
    printf "\nERROR: BlasCh Python script not found at /scripts/submodules/Blasch_PipeCraft.py\n"
    end_process
fi

# Set up BlasCh command arguments
blasch_args="--input_chimeras_dir $input_chimeras_dir"
blasch_args="$blasch_args --self_fasta_dir $self_fasta_dir"
blasch_args="$blasch_args --output_dir $blasch_output_dir"
blasch_args="$blasch_args --threads $threads"

# Add threshold parameters
blasch_args="$blasch_args --high_coverage_threshold $high_coverage_threshold"
blasch_args="$blasch_args --high_identity_threshold $high_identity_threshold"
blasch_args="$blasch_args --borderline_coverage_threshold $borderline_coverage_threshold"
blasch_args="$blasch_args --borderline_identity_threshold $borderline_identity_threshold"

# Add reference database if provided
if [[ "$reference_db" != "undefined" ]] && [[ -n "$reference_db" ]]; then
    if [[ -f "$reference_db" ]]; then
        blasch_args="$blasch_args --reference_db $reference_db"
        printf "Using reference database: $reference_db\n"
    else
        printf "WARNING: reference_db file not found: $reference_db. Proceeding without reference database.\n"
    fi
else
    printf "No reference database specified. Using only self-databases.\n"
fi

printf "\nRunning BlasCh with the following parameters:\n"
printf "Input chimeras directory: $input_chimeras_dir\n"
printf "Self FASTA directory: $self_fasta_dir\n"
printf "Output directory: $blasch_output_dir\n"
printf "Threads: $threads\n"
printf "High coverage threshold: $high_coverage_threshold%%\n"
printf "High identity threshold: $high_identity_threshold%%\n"
printf "Borderline coverage threshold: $borderline_coverage_threshold%%\n"
printf "Borderline identity threshold: $borderline_identity_threshold%%\n"

# Run BlasCh
checkerror=$(python3 /scripts/submodules/Blasch_PipeCraft.py $blasch_args 2>&1)
check_app_error

printf "\nBlasCh analysis completed successfully!\n"
printf "Results are available in: $blasch_output_dir\n"
printf "\nOutput files:\n"
printf "- Rescued sequences: *_non_chimeric.fasta and *_borderline.fasta\n"
printf "- Chimeric sequences: analysis/*_chimeric.fasta\n"
printf "- Multiple alignment sequences: analysis/*_multiple_alignments.fasta\n"
printf "- Analysis report: chimera_detection_report.txt\n"
printf "- Detailed results: analysis/*_sequence_details.csv\n"

#end
printf "\nDone\n"

#!/bin/bash

# Sequence clustering with SWARM incorporating advanced options
# Input = single-end FASTA/FASTQ files.
# Output = Representative OTU sequences (FASTA), SWARM cluster file, and statistics file.

################################################
### Third-party applications:
# swarm (v3.1.5)
# vsearch (for dereplication v2.23.0)
# GNU Parallel (20210422)
# seqkit (v2.3.0)
# pigz
################################################

# Checking tool versions
swarm_version=$(swarm --version 2>&1 | head -n1 | awk '{print $2}')
vsearch_version=$(vsearch --version 2>&1 | head -n1 | awk '{print $2}' | sed 's/,//g')
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# swarm (version %s)\n" "$swarm_version"
printf "# vsearch (version %s)\n" "$vsearch_version"
printf "# seqkit (version %s)\n" "$seqkit_version"

# Load variables (these can be set in the environment or upstream)
# SWARM clustering options with defaults:
swarm_d=${swarm_d:-1}                      # Resolution (differences), default 1.
swarm_no_break=${swarm_no_break}           # If defined (default "true"), adds --no-otu-breaking.
# Fastidious options (only applicable if swarm_d == 1)
swarm_boundary=${swarm_boundary}           # default 3
swarm_ceiling=${swarm_ceiling}             # default 1000 (MB) or leave empty.
swarm_fastidious=${swarm_fastidious}       # Flag, default "true" for -f.
swarm_bloom_bits=${swarm_bloom_bits}       # default 16
# Input/output options:
swarm_append=${swarm_append}               # Value to use when abundance is missing.
swarm_internal_structure=${swarm_internal_structure}  # Filename for internal structure.
swarm_network_file=${swarm_network_file}   # Filename for network output.
swarm_log=${swarm_log}                     # Log filename.
swarm_mothur=${swarm_mothur}               # Flag for mothur-like output.
swarm_uclust_file=${swarm_uclust_file}     # Filename for UCLUST-like output.
# Default output filenames:
swarm_output_file=${swarm_output_file:-amplicon.swarms}
swarm_stats_file=${swarm_stats_file:-amplicon.stats}
swarm_seeds_file=${swarm_seeds_file:-amplicons.fasta}
# Pairwise alignment advanced options (only when swarm_d > 1)
swarm_match=${swarm_match}                 # Default reward for nucleotide match, default 5.
swarm_mismatch=${swarm_mismatch}           # Default penalty for mismatch, default 4.
swarm_gap_open=${swarm_gap_open}           # Gap open penalty, default 12.
swarm_gap_ext=${swarm_gap_ext}             # Gap extension penalty, default 4.
swarm_disable_sse3=${swarm_disable_sse3}   # Flag, default "true" to disable SSE3.
# Thread setting: default to 4 but can be overridden.
swarm_threads=${swarm_threads:-4}

# File format and input file variables
fileFormat=${fileFormat}                   # "fasta", "fastq", etc.
fasta_file=${fasta_file}                   # Single-end input file; if FASTQ, conversion will be performed.

# Source for shared functions
source /scripts/submodules/framework.functions.sh

# Define output directory for swarm clustering
output_dir=$"/input/clustering_out"
export output_dir

###############################
### Multi-run Check
###############################
if [[ -d "/input/multiRunDir" ]]; then
    echo "SWARM clustering pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = clustering"
    cd /input/multiRunDir
    DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "chimeraFiltered_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    echo "Working in directories:"
    echo "$DIRS"
    multiDir="TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = clustering"
    DIRS=$(pwd)
    echo "Working directory: $DIRS"
fi

###############################
### Start of the Workflow
###############################
# Array for storing output FASTA files (for multi-run)
declare -a output_fastas

for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)

    cd "$seqrun"
    if [[ "$multiDir" == "TRUE" ]]; then
        count=$(ls -1 *."$fileFormat" 2>/dev/null | wc -l)
        if [[ $count -eq 0 ]]; then
            printf "[ERROR] No files with extension \"%s\" found in %s. Quitting.\n" "$fileFormat" "$seqrun" >&2
            end_process
        fi
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/clustering_out_swarm"
        export output_dir
        output_fastas+=("$output_dir/swarm_clusters.fasta")
        first_file_check
        prepare_SE_env
    else
        output_dir=$"/input/clustering_out_swarm"
        export output_dir
        output_fastas=("$output_dir/swarm_clusters.fasta")
        first_file_check
        prepare_SE_env
    fi

    ### Pre-process samples: Check files, decompress if needed, and validate formats.
    printf "Checking input files ...\n"
    for file in *."$fileFormat"; do
        check_gz_zip_SE
        check_extension_fastx
    done

    # If input is FASTQ, convert to FASTA
    if [[ "$fileFormat" == "fastq" ]] || [[ "$fileFormat" == "fq" ]]; then
        for file in *."$fileFormat"; do
            samp_name=$(basename "$file" | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
            checkerror=$(seqkit fq2fa -t dna --line-width 0 "$file" -o "$samp_name.fasta" 2>&1)
            check_app_error
        done
        was_fastq="true"
        export was_fastq
        fileFormat="fasta"
        export fileFormat
    fi

    # Create temporary directories for dereplication
    [[ -d tempdir ]] && rm -rf tempdir
    mkdir -p tempdir
    [[ -d dereplicated_sequences ]] && rm -rf dereplicated_sequences
    mkdir -p dereplicated_sequences

    ### Dereplicate individual samples and rename sequences to sha1
    derep_rename () {
        samp_name=$(basename "$1" | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
        vsearch --derep_fulllength "$1" \
                --relabel_sha1 \
                --output - \
                --fasta_width 0 \
                --threads 1 \
                --sizein --sizeout \
        | sed 's/>.*/&;sample='"$samp_name"'/' > dereplicated_sequences/"$samp_name".fasta
    }
    export -f derep_rename
    printf "Dereplicating individual samples ...\n"
    find . -maxdepth 1 -name "*.${fileFormat}" | parallel -j 1 "derep_rename {}"

    ### Global dereplication: Combine individual dereplicated files
    printf "Performing global dereplication ...\n"
    find dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
    | vsearch --derep_fulllength - \
              --output tempdir/Glob_derep.fasta \
              --fasta_width 0 \
              --threads 1 \
              --sizein --sizeout

    ### Build SWARM command options
    swarm_opts="-d ${swarm_d}"
    # Add no-otu-breaking option if defined
    if [[ "$swarm_no_break" == "true" ]]; then
        swarm_opts+=" -n"
    fi
    # If resolution is 1, add fastidious options if provided
    if [[ "$swarm_d" -eq 1 ]]; then
        [[ -n "$swarm_boundary" ]] && swarm_opts+=" -b ${swarm_boundary}"
        [[ -n "$swarm_ceiling" ]] && swarm_opts+=" -c ${swarm_ceiling}"
        if [[ "$swarm_fastidious" == "true" ]]; then
            swarm_opts+=" -f"
        fi
        [[ -n "$swarm_bloom_bits" ]] && swarm_opts+=" -y ${swarm_bloom_bits}"
    else
        # For d > 1, add pairwise alignment advanced options if provided
        [[ -n "$swarm_match" ]] && swarm_opts+=" -m ${swarm_match}"
        [[ -n "$swarm_mismatch" ]] && swarm_opts+=" -p ${swarm_mismatch}"
        [[ -n "$swarm_gap_open" ]] && swarm_opts+=" -g ${swarm_gap_open}"
        [[ -n "$swarm_gap_ext" ]] && swarm_opts+=" -e ${swarm_gap_ext}"
        if [[ "$swarm_disable_sse3" == "true" ]]; then
            swarm_opts+=" -x"
        fi
    fi
    # Append default output options: enable usearch abundance (-z)
    swarm_opts+=" -z"
    swarm_opts+=" -o ${swarm_output_file}"
    swarm_opts+=" -s ${swarm_stats_file}"
    swarm_opts+=" -w ${swarm_seeds_file}"
    # Add threads option (if supported by swarm)
    swarm_opts+=" -t ${swarm_threads}"
    # Additional input/output options if set by user
    [[ -n "$swarm_append" ]] && swarm_opts+=" -a ${swarm_append}"
    [[ -n "$swarm_internal_structure" ]] && swarm_opts+=" -i ${swarm_internal_structure}"
    [[ -n "$swarm_network_file" ]] && swarm_opts+=" -j ${swarm_network_file}"
    [[ -n "$swarm_log" ]] && swarm_opts+=" -l ${swarm_log}"
    [[ "$swarm_mothur" == "true" ]] && swarm_opts+=" -r"
    [[ -n "$swarm_uclust_file" ]] && swarm_opts+=" -u ${swarm_uclust_file}"

    printf "SWARM options: %s\n" "$swarm_opts"

    ### Clustering with SWARM
    printf "Clustering with SWARM ...\n"
    swarm $swarm_opts tempdir/Glob_derep.fasta
    if [[ $? -ne 0 ]]; then
        printf "Error: SWARM clustering failed.\n" >&2
        exit 1
    fi

    #################################################
    ### CLEAN UP AND COMPILE FINAL README FILE
    #################################################
    if [[ "$was_fastq" == "true" ]]; then
        mkdir -p "$output_dir/clustering_input_to_FASTA"
        mv *.fasta "$output_dir/clustering_input_to_FASTA"
    fi
    if [[ $debugger != "true" ]]; then
        rm -rf tempdir dereplicated_sequences
    fi

    ### Create README.txt with clustering details
    count_features "$output_dir/swarm_clusters.fasta"
    end=$(date +%s)
    runtime=$((end - start))
    printf "# Sequences were clustered into swarm_clusters using SWARM clustering with advanced options.
Start time: %s
End time: %s
Runtime: %s seconds

Files in 'clustering_out_swarm' directory:
--------------------------------------------
swarm_clusters.fasta       = FASTA formatted representative OTU sequences (seed sequences).
%s               = SWARM cluster assignment file.
%s               = SWARM statistics file.

Core command ->
swarm %s tempdir/Glob_derep.fasta

NOTE: OTU table generation is not included in this workflow.
" "$start_time" "$(date)" "$runtime" "${swarm_output_file}" "${swarm_stats_file}" "$swarm_opts" > "$output_dir/README.txt"

    printf "\nDONE for run: %s\nTotal time: %s sec.\n" "$seqrun" "$runtime"

    # Return to multiRunDir root if applicable
    if [[ "$multiDir" == "TRUE" ]]; then
        cd /input/multiRunDir
    fi
done

# For multi-run analyses, validate and combine output FASTA files if needed
if [[ "$multiDir" == "TRUE" ]]; then
    valid_fastas=()
    for fasta in "${output_fastas[@]}"; do
        if [[ -f "$fasta" && -s "$fasta" ]]; then
            valid_fastas+=("$fasta")
        else
            printf "Warning: FASTA file not found or empty: %s\n" "$fasta" >&2
        fi
    done
    output_fasta=$(IFS=,; echo "${valid_fastas[*]}")
else
    output_fasta="${output_fastas[0]}"
fi

# Final output variables for downstream services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=single_end"
echo "output_fasta=$output_fasta"
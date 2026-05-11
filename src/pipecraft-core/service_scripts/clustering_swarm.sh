#!/bin/bash

# Sequence clustering with SWARM incorporating advanced options
# Input = single-end FASTA/FASTQ files.
# Output = Representative OTU sequences (FASTA), SWARM cluster file, and statistics file.

################################################
### Third-party applications:
# swarm
# vsearch
# GNU Parallel (20260122)
# seqkit
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
swarm_d=${swarm_d:-1}                     # Resolution (differences), default 1.
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
# Default output filenames (base names only; output_dir is added later):
swarm_output_file=${swarm_output_file:-amplicon.swarms}
swarm_stats_file=${swarm_stats_file:-amplicon.stats}
swarm_seeds_file=${swarm_seeds_file:-swarm_clusters.fasta}
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
output_dir=$"/input/clustering_out_swarm"
export output_dir

###############################
### Single-run mode
###############################
echo "Working with individual sequencing run"
echo "Process = clustering"
DIRS=$(pwd)
echo "Working directory: $DIRS"

###############################
### Start of the Workflow
###############################
# Array for storing output FASTA files (for multi-run)
declare -a output_fastas

for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)

    cd "$seqrun"
    count=$(ls -1 *."$fileFormat" 2>/dev/null | wc -l)
    if [[ $count -eq 0 ]]; then
        printf "[ERROR] No files with extension \"%s\" found in %s. Quitting.\n" "$fileFormat" "$seqrun" >&2
        end_process
    fi
    output_dir=$"/input/clustering_out_swarm"
    export output_dir
    mkdir -p "$output_dir"
    swarm_output_path="$output_dir/$swarm_output_file"
    swarm_stats_path="$output_dir/$swarm_stats_file"
    swarm_seeds_path="$output_dir/$swarm_seeds_file"
    output_fastas=("$swarm_seeds_path")
    first_file_check
    prepare_SE_env

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
    [[ -d tempdir2 ]] && rm -rf tempdir2
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
    # If resolution is 1, add fastidious options ONLY if fastidious mode is enabled
    if [[ "$swarm_d" -eq 1 ]]; then
        if [[ "$swarm_fastidious" == "true" ]]; then
            swarm_opts+=" -f"
            # Add fastidious-specific options only when fastidious mode is enabled
            [[ -n "$swarm_boundary" ]] && swarm_opts+=" -b ${swarm_boundary}"
            [[ -n "$swarm_ceiling" ]] && swarm_opts+=" -c ${swarm_ceiling}"
            [[ -n "$swarm_bloom_bits" ]] && swarm_opts+=" -y ${swarm_bloom_bits}"
        fi
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
    swarm_opts+=" -o ${swarm_output_path}"
    swarm_opts+=" -s ${swarm_stats_path}"
    swarm_opts+=" -w ${swarm_seeds_path}"
    # Add threads option
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
    ### OTU Table Generation from SWARM output (Polars/NumPy-accelerated)
    printf "Generating OTU table ...\n"
    otu_table_path="$output_dir/swarm_clusters_table.txt"

    python3 - "dereplicated_sequences" \
              "$swarm_output_path" \
              "$swarm_seeds_path" \
              "$otu_table_path" <<'PYEOF'
import sys, os, re
import numpy as np

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

derep_dir   = sys.argv[1]
swarms_file = sys.argv[2]
seeds_file  = sys.argv[3]
out_file    = sys.argv[4]

# -------------------------------------------------------
# STEP 1: Bulk-parse all per-sample derep FASTAs at once.
# A single sha1 can occur in multiple samples (metabarcoding); we must
# record EVERY (sample, size) pair so the OTU table credits all contributors.
# -------------------------------------------------------
HEADER_RE = re.compile(
    rb'>([0-9a-f]{40});size=(\d+);sample=([^;\s\n]+)',
    re.MULTILINE
)

amp_to_samples = {}   # sha1(str) -> list[(sample(str), size(int))]
sample_set = set()

for fname in sorted(os.listdir(derep_dir)):
    if not fname.endswith('.fasta'):
        continue
    fpath = os.path.join(derep_dir, fname)
    with open(fpath, 'rb') as fh:
        blob = fh.read()
    for m in HEADER_RE.finditer(blob):
        sha    = m.group(1).decode()
        size   = int(m.group(2))
        sample = m.group(3).decode()
        amp_to_samples.setdefault(sha, []).append((sample, size))
        sample_set.add(sample)

# -------------------------------------------------------
# STEP 2: Sample index
# -------------------------------------------------------
samples    = sorted(sample_set)
sample_idx = {s: i for i, s in enumerate(samples)}
n_samples  = len(samples)
print(f"Samples: {n_samples}  Amplicons: {len(amp_to_samples)}", flush=True)

# -------------------------------------------------------
# STEP 3: Seed IDs in order from seeds FASTA
# -------------------------------------------------------
seed_ids = []
SEED_RE  = re.compile(rb'^>([^;\s]+)', re.MULTILINE)
with open(seeds_file, 'rb') as fh:
    blob = fh.read()
for m in SEED_RE.finditer(blob):
    seed_ids.append(m.group(1).decode())

# -------------------------------------------------------
# STEP 4: Parse swarms file + accumulate counts into numpy
# Byte-level parsing to avoid UTF-8 decode overhead per entry
# -------------------------------------------------------
otu_ids   = []
count_mat = []

with open(swarms_file, 'rb') as fh:
    for otu_idx, line in enumerate(fh):
        line = line.rstrip()
        if not line:
            continue

        counts = np.zeros(n_samples, dtype=np.int32)

        for amp_entry in line.split(b' '):
            sha = amp_entry.split(b';', 1)[0].decode()
            for samp, sz in amp_to_samples.get(sha, ()):
                counts[sample_idx[samp]] += sz

        otu_id = seed_ids[otu_idx] if otu_idx < len(seed_ids) else f"OTU_{otu_idx}"
        otu_ids.append(otu_id)
        count_mat.append(counts)

# Stack all OTU rows into a single 2D C-contiguous numpy array
count_arr = np.vstack(count_mat)
n_otus    = len(otu_ids)
print(f"OTUs: {n_otus}", flush=True)

# -------------------------------------------------------
# STEP 5: Write TXT OTU table with header: "Swarm_cluster_ID" + sample names
# Polars path: build DataFrame from numpy arrays (zero-copy via Arrow)
# Fallback: write with python buffered I/O
# -------------------------------------------------------
if HAS_POLARS:
    col_dict = {"Swarm_cluster_ID": otu_ids}
    for i, s in enumerate(samples):
        col_dict[s] = count_arr[:, i]   # numpy column slice, zero-copy into Arrow
    df = pl.DataFrame(col_dict)
    df.write_csv(out_file, separator="\t")
    print(f"Written via Polars (zero-copy Arrow build)", flush=True)
else:
    # Fallback: buffered write with pre-joined lines
    header = "Swarm_cluster_ID\t" + "\t".join(samples) + "\n"
    lines  = [header]
    for i, otu_id in enumerate(otu_ids):
        lines.append(otu_id + "\t" + "\t".join(map(str, count_arr[i])) + "\n")
        if len(lines) >= 50000:             # flush in chunks to avoid RAM spike
            with open(out_file, 'a') as fh:
                fh.writelines(lines)
            lines = []
    if lines:
        mode = 'w' if not os.path.exists(out_file) else 'a'
        with open(out_file, mode) as fh:
            fh.writelines(lines)

print(f"Done: {n_otus} OTUs x {n_samples} samples -> {out_file}", flush=True)
PYEOF

    if [[ $? -ne 0 ]]; then
        printf "Error: OTU table generation failed.\n" >&2
        exit 1
    fi
    printf "OTU table written to: %s\n" "$otu_table_path"

    #################################################
    ### CLEAN UP AND COMPILE FINAL README FILE
    #################################################
    ### Create README.txt with clustering details
    input_file_count=$(ls -1 *."$fileFormat" 2>/dev/null | wc -l)
    if [[ $input_file_count -gt 0 ]]; then
        total_reads=$(seqkit stats -T *."$fileFormat" | awk 'NR>1 {sum+=$4} END {print sum+0}')
    else
        total_reads=0
    fi
    if [[ "$was_fastq" == "true" ]]; then
        mkdir -p "$output_dir/clustering_input_to_FASTA"
        mv *.fasta "$output_dir/clustering_input_to_FASTA"
    fi
    if [[ $debugger != "true" ]]; then
        rm -rf tempdir tempdir2 dereplicated_sequences
    fi
    cluster_count=$(grep -c "^>" "$swarm_seeds_path" 2>/dev/null || echo 0)
    end=$(date +%s)
    runtime=$((end - start))
    seeds_basename=$(basename "$swarm_seeds_path")
    swarms_basename=$(basename "$swarm_output_path")
    stats_basename=$(basename "$swarm_stats_path")
    core_command="$swarm_opts"
    if [[ "$output_dir" == /input/* ]]; then
        relative_output_dir="${output_dir#/input/}"
        core_command="${core_command//$output_dir/$relative_output_dir}"
    fi
    printf "# Reads were clustered using SWARM (see 'Core command' below for the used settings).
Start time: %s
End time: %s
Runtime: %s seconds

Files in 'clustering_out_swarm' directory:
--------------------------------------------
%s = FASTA formatted representative swarm-cluster sequences (seed sequences).
%s = SWARM cluster assignment file.
%s = SWARM statistics file.
swarm_clusters_table.txt = SWARM clusters table (tab delimited file).

Input files processed               = %s
Total input reads                   = %s
Number of seed sequences (clusters) = %s

Core command ->
swarm %s tempdir/Glob_derep.fasta

" "$start_time" "$(date)" "$runtime" "$seeds_basename" "$swarms_basename" "$stats_basename" "$input_file_count" "$total_reads" "$cluster_count" "$core_command" > "$output_dir/README.txt"

    if [[ "$was_fastq" == "true" ]]; then
        printf "\n\nInput was fastq; converted those to fasta before clustering.\nConverted fasta files in directory 'clustering_input_to_FASTA'\n" >> "$output_dir/README.txt"
    fi

    printf "\n\n##############################################\n###Third-party applications for this process:\n#swarm (version $swarm_version)\n    #citation: Swarm v3: towards tera-scale amplicon clustering. Mahé F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T. (2021) Bioinformatics doi: 10.1093/bioinformatics/btab493\n    #https://github.com/torognes/swarm\n#vsearch (version $vsearch_version)\n    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584\n    #https://github.com/torognes/vsearch\n#GNU Parallel (version 20260122)\n    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607\n################################################\n" >> "$output_dir/README.txt"

    printf "\nDONE for run: %s\nTotal time: %s sec.\n" "$seqrun" "$runtime"

    # end single-run iteration
done

# Single-run: set output_fasta to the generated seeds file
output_fasta="${output_fastas[0]}"

# Final output variables for downstream services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=single_end"
echo "output_fasta=$output_fasta"
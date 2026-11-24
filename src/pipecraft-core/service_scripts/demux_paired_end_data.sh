#!/bin/bash

# Demultiplex PAIRED-END reads
# Demultiplexing of paired-end reads in mixed orientation using single-end or paired-end indexes is supported.
# Input = a directory with fastq/fasta files (R1.fastq; R2.fastq); and indexes file in fasta format (header as a sample name).
#
# Dual indexes are demultiplexed in two steps so cutadapt only opens files for
# index combinations that are listed in the indexes file:
#   1) unique F indexes
#   2) for each F bin, only the R indexes paired with that F
# Mixed orientation is still handled as two cutadapt rounds (R1/R2 swapped).

################################################
###Third-party applications:
#cutadapt
#seqkit
#pigz
#python3 with biopython
##############################################
#Checking tool versions
cutadapt_version=$(cutadapt --version 2>&1)
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# cutadapt (version $cutadapt_version)\n"
printf "# seqkit (version $seqkit_version)\n"

#Load variables
regex='[^/]*$'
oligos_file_path=$(echo $index_file | grep -oP "$regex")
oligos_file=$(basename $oligos_file_path) #basename, needed for macOS
indexes_file=$(printf "/extraFiles/$oligos_file")
error_rate="-e ${index_mismatch}"
if [ "$no_indels" = true ] ; then
    no_indels=$"--no-indels"
else
    no_indels=''
fi
minlen=$"--minimum-length ${min_length}"
overlap=$"--overlap ${overlap}"

### Check CPU cores ###
# 'cores' is passed in by the app (Resource Manager CPU setting). Validate it,
# then compare against the cores actually available inside the container
# (nproc), using the smaller valid value so we never request more cores than
# exist or run with a missing/garbage value. This single 'cores' value is then
# used for every cutadapt (--cores) and seqkit (--threads) call below.
detected_cores=$(nproc 2>/dev/null || echo 1)
if [[ "$cores" =~ ^[0-9]+$ ]] && (( cores >= 1 )); then
    if (( cores > detected_cores )); then
        printf "# WARNING: requested %s cores but container has %s; using %s\n" "$cores" "$detected_cores" "$detected_cores"
        cores=$detected_cores
    fi
else
    printf "# WARNING: invalid 'cores' value ('%s'); using detected %s\n" "$cores" "$detected_cores"
    cores=$detected_cores
fi
printf "# Using %s CPU core(s) for cutadapt and seqkit\n" "$cores"

printf "overlap = $overlap\n"
printf "no_indels = $no_indels\n"
printf "error_rate = $error_rate\n"
printf "search_window = $search_window\n"
printf "minlen = $minlen\n"
printf "cores = $cores\n"

# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/demultiplex_out"
#python module: per unique F, write only the R indexes listed for that F
build_allowed_R_per_F=$"python3 /scripts/submodules/build_allowed_R_per_F.demuxModule.py"

# Soft-raise the open-files limit. Rootless Podman often cannot raise it;
# cutadapt then uses whatever the user namespace already allows.
ulimit -S -n 6000 2>/dev/null || true

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env #checks also multiple R1/R2 occurrences
### Check barcodes file
check_indexes_file

cutadapt_demux () {
    checkerror=$(cutadapt --quiet \
        $error_rate \
        $no_indels \
        $overlap \
        $minlen \
        --cores ${cores} \
        "$@" 2>&1)
    check_app_error
}

prepare_dual_index_files () {
    printf "Preparing dual-index files for two-step demux ...\n"
    sed -e 's/\.\.\..*//' < tempdir2/ValidatedBarcodesFileForDemux.fasta.temp > tempdir2/index_fwd.fasta

    # Uppercase before rmdup so the same F sequence is not split by case
    # (Python matches F sequences case-insensitively).
    seqkit seq --quiet -u -w 0 tempdir2/index_fwd.fasta | \
        seqkit rmdup --quiet --by-seq -w 0 > tempdir2/index_fwd.uniq.fasta
    seqkit replace --quiet tempdir2/index_fwd.uniq.fasta -w 0 -p .+ -r "indexF_{nr}" > tempdir2/index_fwd.uniq.renamed.fasta

    sed -e '/^>/!s/^/search_window/' tempdir2/index_fwd.uniq.renamed.fasta > $output_dir/index_fwd.fasta
    sed -i "s/search_window/XN{$search_window}/" $output_dir/index_fwd.fasta

    mkdir -p tempdir2/R_per_F
    $build_allowed_R_per_F \
        tempdir2/ValidatedBarcodesFileForDemux.fasta.temp \
        tempdir2/index_fwd.uniq.renamed.fasta \
        tempdir2/R_per_F \
        "$search_window"
    if [[ $? -ne 0 ]]; then
        printf '%s\n' "ERROR]: failed to build per-F reverse-index files from the indexes file.
>Quitting" >&2
        end_process
    fi
    if ! compgen -G "tempdir2/R_per_F/R_for_*.fasta" > /dev/null; then
        printf '%s\n' "ERROR]: no per-F reverse-index files were written.
>Quitting" >&2
        end_process
    fi
    cat tempdir2/R_per_F/R_for_*.fasta > $output_dir/index_rev.fasta
}

prepare_single_index_file () {
    sed -i '/^>/!s/^/search_window/' tempdir2/ValidatedBarcodesFileForDemux.fasta.temp
    sed -i "s/search_window/XN{$search_window}/" tempdir2/ValidatedBarcodesFileForDemux.fasta.temp
    mv tempdir2/ValidatedBarcodesFileForDemux.fasta.temp $output_dir/index_file.fasta
}

merge_round_pair () {
    local round1_file=$1
    local round2_file=$2
    local final_file=$3
    if [[ -f "$round1_file" && -f "$round2_file" ]]; then
        cat "$round1_file" "$round2_file" > "$final_file"
        rm -f "$round1_file" "$round2_file"
    elif [[ -f "$round1_file" ]]; then
        mv "$round1_file" "$final_file"
    elif [[ -f "$round2_file" ]]; then
        mv "$round2_file" "$final_file"
    fi
}

merge_demux_rounds () {
    # Merge round1 + round2 per sample (a sample may be present in only one round)
    local dir=$1
    local samples_file=$dir/demux_samples_to_merge.txt
    local f base sample
    : > "$samples_file"
    for f in "$dir"/round1-*.R1."$fileFormat" "$dir"/round2-*.R1."$fileFormat"; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        if [[ "$base" == round1-* ]]; then
            sample=${base#round1-}
        else
            sample=${base#round2-}
        fi
        sample=${sample%.R1.$fileFormat}
        if [[ "$sample" == "unknown" || "$sample" == "unknown-unknown" ]]; then
            continue
        fi
        printf "%s\n" "$sample" >> "$samples_file"
    done
    if [[ -s "$samples_file" ]]; then
        while read -r sample; do
            [[ -n "$sample" ]] || continue
            merge_round_pair \
                "$dir/round1-${sample}.R1.$fileFormat" \
                "$dir/round2-${sample}.R1.$fileFormat" \
                "$dir/${sample}.R1.$fileFormat"
            merge_round_pair \
                "$dir/round1-${sample}.R2.$fileFormat" \
                "$dir/round2-${sample}.R2.$fileFormat" \
                "$dir/${sample}.R2.$fileFormat"
        done < <(sort -u "$samples_file")
    fi
    rm -f "$samples_file"
}

append_to_unknown () {
    local src_r1=$1
    local src_r2=$2
    if [[ ! -s "$src_r1" ]]; then
        rm -f "$src_r1" "$src_r2"
        return 0
    fi
    cat "$src_r1" >> $output_dir/unknown.R1.$fileFormat
    cat "$src_r2" >> $output_dir/unknown.R2.$fileFormat
    rm -f "$src_r1" "$src_r2"
}

assign_allowed_R_per_F () {
    # $1 = round prefix (round1|round2)
    # $2 = orientation: normal (R on R2, swap I/O) or rc (R on R1)
    local round_prefix=$1
    local orientation=$2
    local rfa fname r1 r2
    for rfa in tempdir2/R_per_F/R_for_*.fasta; do
        [[ -f "$rfa" ]] || continue
        fname=$(basename "$rfa" .fasta)
        fname=${fname#R_for_}
        r1=$output_dir/byF/${round_prefix}-${fname}.R1.$fileFormat
        r2=$output_dir/byF/${round_prefix}-${fname}.R2.$fileFormat
        if [[ ! -s "$r1" || ! -s "$r2" ]]; then
            continue
        fi
        printf "   %s %s: assigning samples by allowed R indexes ...\n" "$round_prefix" "$fname"
        if [[ "$orientation" == "normal" ]]; then
            # R index is on R2; demux keys off the first input file, so swap
            cutadapt_demux \
                -g "file:$rfa" \
                -o $output_dir/${round_prefix}-{name}.R2.$fileFormat \
                -p $output_dir/${round_prefix}-{name}.R1.$fileFormat \
                --untrimmed-output $output_dir/unnamed/${round_prefix}-${fname}-unassigned.R2.$fileFormat \
                --untrimmed-paired-output $output_dir/unnamed/${round_prefix}-${fname}-unassigned.R1.$fileFormat \
                $r2 $r1
        else
            # RC round: R index is on original R1 (already the .R1 file)
            cutadapt_demux \
                -g "file:$rfa" \
                -o $output_dir/${round_prefix}-{name}.R1.$fileFormat \
                -p $output_dir/${round_prefix}-{name}.R2.$fileFormat \
                --untrimmed-output $output_dir/unnamed/${round_prefix}-${fname}-unassigned.R1.$fileFormat \
                --untrimmed-paired-output $output_dir/unnamed/${round_prefix}-${fname}-unassigned.R2.$fileFormat \
                $r1 $r2
        fi
        rm -f "$r1" "$r2"
    done
}

demux_dual_two_step () {
    mkdir -p $output_dir/byF $output_dir/unnamed

    printf "   Round1: demultiplex by unique F indexes ...\n"
    cutadapt_demux \
        -g file:$output_dir/index_fwd.fasta \
        -o $output_dir/byF/round1-{name}.R1.$fileFormat \
        -p $output_dir/byF/round1-{name}.R2.$fileFormat \
        $inputR1.$fileFormat $inputR2.$fileFormat

    assign_allowed_R_per_F "round1" "normal"

    if [[ -s $output_dir/byF/round1-unknown.R1.$fileFormat && -s $output_dir/byF/round1-unknown.R2.$fileFormat ]]; then
        printf "   Round2 (RC; R1 and R2 position switched): demultiplex by unique F indexes ...\n"
        cutadapt_demux \
            -g file:$output_dir/index_fwd.fasta \
            -o $output_dir/byF/round2-{name}.R2.$fileFormat \
            -p $output_dir/byF/round2-{name}.R1.$fileFormat \
            $output_dir/byF/round1-unknown.R2.$fileFormat \
            $output_dir/byF/round1-unknown.R1.$fileFormat
        rm -f $output_dir/byF/round1-unknown.R1.$fileFormat $output_dir/byF/round1-unknown.R2.$fileFormat
        assign_allowed_R_per_F "round2" "rc"
        append_to_unknown \
            $output_dir/byF/round2-unknown.R1.$fileFormat \
            $output_dir/byF/round2-unknown.R2.$fileFormat
    else
        append_to_unknown \
            $output_dir/byF/round1-unknown.R1.$fileFormat \
            $output_dir/byF/round1-unknown.R2.$fileFormat
    fi

    merge_demux_rounds "$output_dir"

    # F matched but R not in the allowed list for that F -> unknown.R1/R2
    local unassigned_r1 unassigned_r2
    for unassigned_r1 in "$output_dir"/unnamed/*-unassigned.R1."$fileFormat"; do
        [[ -f "$unassigned_r1" ]] || continue
        unassigned_r2=${unassigned_r1%.R1.$fileFormat}.R2.$fileFormat
        append_to_unknown "$unassigned_r1" "$unassigned_r2"
    done

    # F bins that were never R-assigned (empty or skipped) -> unknown, not deleted
    local leftover_r1 leftover_r2
    for leftover_r1 in "$output_dir"/byF/*.R1."$fileFormat"; do
        [[ -f "$leftover_r1" ]] || continue
        leftover_r2=${leftover_r1%.R1.$fileFormat}.R2.$fileFormat
        append_to_unknown "$leftover_r1" "$leftover_r2"
    done
    rm -rf $output_dir/byF $output_dir/unnamed
}

demux_single_two_round () {
    printf "   Round1: search index on R1 ...\n"
    cutadapt_demux \
        -g file:$output_dir/index_file.fasta \
        -o $output_dir/round1-{name}.R1.$fileFormat \
        -p $output_dir/round1-{name}.R2.$fileFormat \
        $inputR1.$fileFormat $inputR2.$fileFormat

    if [[ -s $output_dir/round1-unknown.R1.$fileFormat && -s $output_dir/round1-unknown.R2.$fileFormat ]]; then
        printf "   Round2 (RC; R1 and R2 position switched): search index on leftover R2 ...\n"
        cutadapt_demux \
            -g file:$output_dir/index_file.fasta \
            -o $output_dir/round2-{name}.R2.$fileFormat \
            -p $output_dir/round2-{name}.R1.$fileFormat \
            $output_dir/round1-unknown.R2.$fileFormat \
            $output_dir/round1-unknown.R1.$fileFormat
        rm -f $output_dir/round1-unknown.R1.$fileFormat $output_dir/round1-unknown.R2.$fileFormat
        append_to_unknown \
            $output_dir/round2-unknown.R1.$fileFormat \
            $output_dir/round2-unknown.R2.$fileFormat
    else
        append_to_unknown \
            $output_dir/round1-unknown.R1.$fileFormat \
            $output_dir/round1-unknown.R2.$fileFormat
    fi

    merge_demux_rounds "$output_dir"
}

if [[ $tag == "dual" ]]; then
    prepare_dual_index_files
else
    prepare_single_index_file
fi

### Process file
printf "Checking files ...\n"
while read LINE; do
    #Write file name without extension
    inputR1=$(echo $LINE | sed -e "s/.$fileFormat//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')

    printf "\n# Demultiplexing with $tag indexes ... \n"
    if [[ $tag == "dual" ]]; then
        demux_dual_two_step
    else
        demux_single_two_round
    fi
done < tempdir2/paired_end_files.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats_demux

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
readme_footer="
Summary of sequence counts in 'seq_count_summary.txt'

################################################
###Third-party applications used for this process:
#cutadapt (version $cutadapt_version) for demultiplexing
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit (version $seqkit_version) for validating indexes file and adjusting sample names
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##############################################"

if [[ $tag == "dual" ]]; then
    printf "# Demultiplexing was performed using cutadapt (paired-end / dual indexes; see 'Core command' below).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Indexes file: $oligos_file (paired-end indexes, FWD...REV per sample).
index_fwd.fasta / index_rev.fasta = unique F indexes and listed R indexes with search window size for cutadapt.

Paired-end data were demultiplexed in two steps so only index combinations listed in the indexes file are written:
  1) unique F indexes
  2) for each F bin, only the R indexes paired with that F
Mixed orientation was handled with two cutadapt rounds (R1/R2 swapped in round 2).
Output R1 and R2 reads are synchronized for merging paired-end data.

Reads that could not be assigned to a listed index combination (no matching F index,
or F matched but the R index is not listed for that F) are in unknown.R1/R2.

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.

Core commands ->
Round1 F: cutadapt -g file:index_fwd.fasta $error_rate $no_indels $overlap $minlen -o byF/round1-{name}.R1 -p byF/round1-{name}.R2 inputR1 inputR2
Round1 R (per unique F; R1/R2 swapped so the R index is on the first file): cutadapt -g file:R_for_{F}.fasta $error_rate $no_indels $overlap $minlen -o round1-{name}.R2 -p round1-{name}.R1 byF/round1-{F}.R2 byF/round1-{F}.R1
Round2 F (RC; R1 and R2 position switched!): cutadapt -g file:index_fwd.fasta $error_rate $no_indels $overlap $minlen -o byF/round2-{name}.R2 -p byF/round2-{name}.R1 round1-unknown.R2 round1-unknown.R1
Round2 R (per unique F): cutadapt -g file:R_for_{F}.fasta $error_rate $no_indels $overlap $minlen -o round2-{name}.R1 -p round2-{name}.R2 byF/round2-{F}.R1 byF/round2-{F}.R2
" > $output_dir/README.txt
else
    printf "# Demultiplexing was performed using cutadapt (single-end indexes; see 'Core command' below).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Indexes file: $oligos_file (single-end indexes; one barcode sequence per sample).
index_file.fasta = $oligos_file with added search window size for cutadapt.

Paired-end data were demultiplexed by searching the index on R1 (round 1) and, for leftover reads, on R2 with R1/R2 swapped (round 2; reverse-complement orientation).
Output R1 and R2 reads are synchronized for merging paired-end data.

Reads with no matching index are in unknown.R1/R2.

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.

Core commands ->
Round1: cutadapt -g file:index_file.fasta $error_rate $no_indels $overlap $minlen -o round1-{name}.R1 -p round1-{name}.R2 inputR1 inputR2
Round2 (RC; R1 and R2 position switched!): cutadapt -g file:index_file.fasta $error_rate $no_indels $overlap $minlen -o round2-{name}.R2 -p round2-{name}.R1 round1-unknown.R2 round1-unknown.R1
" > $output_dir/README.txt
fi

printf "%s" "$readme_footer" >> $output_dir/README.txt

###Done, files in $output_dir folder
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"

#!/bin/bash

# Demultiplex SINGLE-END reads.
# Demultiplexing of single-end reads in mixed orientation using paired-end or single-end indexes is supported.
# Input = a directory with fastq/fasta files; and indexes file in fasta format (header as a sample name).

################################################
###Third-party applications:
#cutadapt
#seqkit
#pigz
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
###############################
###############################

# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/demultiplex_out"

prepare_dual_linked_index_file () {
    printf "Preparing dual-index linked adapters (FWD 5', RC(REV) 3') ...\n"
    sed -e 's/\.\.\..*//' < tempdir2/ValidatedBarcodesFileForDemux.fasta.temp > tempdir2/index_fwd.fasta

    seqkit seq --quiet -n tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
        sed -e 's/^/>/' > tempdir2/sample_names.txt
    awk '!/^>/ {
        p = index($0, "...")
        if (p == 0) {
            print "ERROR]: dual-index sequence is missing FWD...REV" > "/dev/stderr"
            exit 1
        }
        print substr($0, p + 3)
    }' tempdir2/ValidatedBarcodesFileForDemux.fasta.temp > tempdir2/index_rev.temp
    if [[ $? -ne 0 ]]; then
        printf '%s\n' "ERROR]: failed to extract reverse indexes from the indexes file.
>Quitting" >&2
        end_process
    fi

    n_names=$(grep -c '^>' tempdir2/sample_names.txt)
    n_rev=$(grep -cve '^$' tempdir2/index_rev.temp)
    if [[ "$n_names" -ne "$n_rev" ]]; then
        printf '%s\n' "ERROR]: number of sample names ($n_names) does not match reverse indexes ($n_rev).
>Quitting" >&2
        end_process
    fi
    paste -d '\n' tempdir2/sample_names.txt tempdir2/index_rev.temp > tempdir2/index_rev.fasta
    rm -f tempdir2/index_rev.temp tempdir2/sample_names.txt

    checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/index_rev.fasta > tempdir2/index_revRC.fasta 2>&1)
    check_app_error

    # Linked adapters: 5' FWD with start-window, 3' RC(REV) with end-window
    tr "\n" "\t" < tempdir2/index_fwd.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > tempdir2/index_fwd.temp
    tr "\n" "\t" < tempdir2/index_revRC.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > tempdir2/index_revRC.temp
    sed -i "s/\t/\tXN{$search_window}/" tempdir2/index_fwd.temp
    sed -i "s/\t$/XN{$search_window}/" tempdir2/index_revRC.temp
    awk 'BEGIN {FS=OFS="\t"} FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2}' \
        tempdir2/index_fwd.temp tempdir2/index_revRC.temp > tempdir2/paired_index.temp
    sed -e 's/\t/\n/' < tempdir2/paired_index.temp | sed -e 's/\t/\.\.\./' > tempdir2/index_file.fasta
    mv tempdir2/index_file.fasta $output_dir
}

prepare_single_index_file () {
    sed -i '/^>/!s/^/search_window/' tempdir2/ValidatedBarcodesFileForDemux.fasta.temp
    sed -i "s/search_window/XN{$search_window}/" tempdir2/ValidatedBarcodesFileForDemux.fasta.temp
    mv tempdir2/ValidatedBarcodesFileForDemux.fasta.temp $output_dir/index_file.fasta
}

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
### Check barcodes file
check_indexes_file

if [[ $tag == "dual" ]]; then
    prepare_dual_linked_index_file
else
    prepare_single_index_file
fi
indexes_file_in=$"-g file:$output_dir/index_file.fasta"
out=$"-o $output_dir/{name}.$fileFormat"

printf "Checking the input file ...\n"
file_count=$(grep -cve '^$' tempdir2/files_in_folder.txt)
if (( file_count > 1 )); then
    printf '%s\n' "ERROR]: please include only one $fileFormat file in the WORKDIR
>Quitting" >&2
    end_process
fi

while read -r file; do
    [[ -n "$file" && -f "$file" ]] || continue
    input="${file%.$fileFormat}"

    printf "\n# Demultiplexing with $tag indexes ... \n"
    checkerror=$(cutadapt --quiet \
        $indexes_file_in \
        $error_rate \
        $no_indels \
        --revcomp \
        --untrimmed-output $output_dir/unknown.$fileFormat \
        $overlap \
        $minlen \
        --cores ${cores} \
        $out \
        $input.$fileFormat 2>&1)
    check_app_error
done < tempdir2/files_in_folder.txt

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

##############################################
###Third-party applications for this process:
#cutadapt (version $cutadapt_version) for demultiplexing
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit (version $seqkit_version) for validating indexes file
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##############################################"

if [[ $tag == "dual" ]]; then
    printf "# Demultiplexing was performed using cutadapt (single-end reads / dual indexes; see 'Core command' below).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Indexes file: $oligos_file (paired-end indexes, FWD...REV per sample).
index_file.fasta = linked adapters per listed sample: XN{window}FWD...RC(REV)XN{window}
  (5' FWD in the start window; 3' reverse-complemented REV in the end window).

Mixed orientation was handled with --revcomp. Reads that matched on the reverse complement
were reverse complemented (sequence name appended with 'rc'), so output is 5'-3'
relative to FWD...REV. Unassigned reads are in unknown.

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.

Core command ->
cutadapt -g file:index_file.fasta $error_rate $no_indels --revcomp --untrimmed-output unknown $overlap $minlen -o {name} input
" > $output_dir/README.txt
else
    printf "# Demultiplexing was performed using cutadapt (single-end reads / single-end indexes; see 'Core command' below).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Indexes file: $oligos_file (single-end indexes; one barcode sequence per sample).
index_file.fasta = $oligos_file with XN{window} added at the 5' end of each index.

Mixed orientation was handled with --revcomp. Reads that matched on the reverse complement
were reverse complemented (sequence name appended with 'rc'), so all sequences are in
uniform orientation. Unassigned reads are in unknown.

Sequence orientation in demultiplex_out:
  - index on the 5' end of the original read -> output is 5'-3'
  - index on the 3' end of the original read -> the read is reverse complemented; remaining sequence is 3'-5' relative to the original molecule

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.

Core command ->
cutadapt -g file:index_file.fasta $error_rate $no_indels --revcomp --untrimmed-output unknown $overlap $minlen -o {name} input
" > $output_dir/README.txt
fi

printf "%s" "$readme_footer" >> $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=single_end"

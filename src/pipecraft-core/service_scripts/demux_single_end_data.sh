#!/bin/bash

# Demultiplex SINGLE-END reads.
# Demultiplexing of single-end reads in mixed orientation using paired-end or single-end indexes is supported.
# Input = a directory with fastq/fasta files; and indexes file in fasta format (header as a sample name).

##########################################################
###Third-party applications:
#cutadapt v4.4
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #Distributed under the MIT license
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz v2.4
##################################################################

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
minlen=$"--minimum-length ${min_seq_length}"
cores=$"--cores ${cores}"
overlap=$"--overlap ${overlap}"
search_window=${search_window}
###############################
###############################

# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/demultiplex_out"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
### Check barcodes file
check_indexes_file

### Process file
printf "Checking the input file ...\n"
#Chech that only one $fileFormat file is in the WORRKING dir
files=$(ls $workingDir | grep -c ".$fileFormat")
if (( $files > 1 )); then
    printf '%s\n' "ERROR]: please include only one $fileFormat file in the WORKDIR \n
>Quitting" >&2
    end_process
fi

for file in *.$fileFormat; do
    #Write file name without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_gz_zip_SE
    ### Check input formats (fastq/fasta supported)
    check_extension_fastx

    ### Check if dual indexes or single indexes and prepare workflow accordingly
    if grep -q "\..." tempdir2/ValidatedBarcodesFileForDemux.fasta.temp; then
        #dual indexes
        #make rev indexes file
        seqkit seq --quiet -n tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
        sed -e 's/^/>/' > tempdir2/sample_names.txt
        grep "\..." tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
        awk 'BEGIN{FS="."}{print $4}' > tempdir2/index_rev.temp
        touch tempdir2/index_rev.fasta
        i=1
        p=$"p"
        while read HEADER; do
            echo $HEADER >> tempdir2/index_rev.fasta
            sed --quiet $i$p tempdir2/index_rev.temp >> tempdir2/index_rev.fasta
            i=$(($i + 1))
        done < tempdir2/sample_names.txt
        rm tempdir2/index_rev.temp
        #make fwd indexes file
        sed -e 's/\.\.\..*//' < tempdir2/ValidatedBarcodesFileForDemux.fasta.temp > tempdir2/index_fwd.fasta
        #reverse complementary REV indexes
        checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/index_rev.fasta > tempdir2/index_revRC.fasta 2>&1)
        check_app_error
        #Make linked indexes files where REV indexes are in RC orientation
        tr "\n" "\t" < tempdir2/index_fwd.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > tempdir2/index_fwd.temp
        tr "\n" "\t" < tempdir2/index_revRC.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > tempdir2/index_revRC.temp
        sed -i "s/\t/\tXN{$search_window}/" tempdir2/index_fwd.temp #add search window size to indexes
        sed -i "s/\t$/XN{$search_window}/" tempdir2/index_revRC.temp #add search window size to indexes
        awk 'BEGIN {FS=OFS="\t"} FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2}' tempdir2/index_fwd.temp tempdir2/index_revRC.temp > tempdir2/paired_index.temp
        sed -e 's/\t/\n/' < tempdir2/paired_index.temp | sed -e 's/\t/\.\.\./' > tempdir2/index_file.fasta

        #assign demux variables
        #REV indexes are 5'-3' orientation for cutadapt search 
        mv tempdir2/index_file.fasta $output_dir #move edited indexes file with window size into output_dir 
        indexes_file_in=$"-g file:$output_dir/index_file.fasta"
        out=$"-o $output_dir/{name}.$extension"
    else
        #single indexes
        # Add search window size to indexes 
        sed -i '/^>/!s/^/search_window/' tempdir2/ValidatedBarcodesFileForDemux.fasta.temp 
        sed -i "s/search_window/XN{$search_window}/" tempdir2/ValidatedBarcodesFileForDemux.fasta.temp 
        #Move edited indexes file with window size into output_dir 
        mv tempdir2/ValidatedBarcodesFileForDemux.fasta.temp tempdir2/index_file.fasta
        mv tempdir2/index_file.fasta $output_dir
        #assign demux variables
        indexes_file_in=$"-g file:$output_dir/index_file.fasta" 
        out=$"-o $output_dir/{name}.$extension"
    fi

    ############################
    ### Start demultiplexing ###
    ############################
    printf "\n# Demultiplexing with $tag indexes ... \n"
    ### Demultiplex with cutadapt
    checkerror=$(cutadapt --quiet \
    $indexes_file_in \
    $error_rate \
    $no_indels \
    --revcomp \
    --untrimmed-output $output_dir/unknown.$extension \
    $overlap \
    $minlen \
    $cores \
    $out \
    $input.$extension 2>&1)
    check_app_error
done

#Remove 'rc' string from seq if the indexes were found on reverse complementary strand
for file in $output_dir/*.$extension; do
    awk -i inplace '{ gsub(/ rc$/, "") }; { print }' $file
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats_demux
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
printf "# Demultiplexing was performed using cutadapt (see 'Core command' below for the used settings).

Files in 'demultiplex_out' directory represent per sample sequence files, that were generated based on the specified indexes file ($oligos_file).
index_file.fasta = $oligos_file but with added search window size for cutadapt.

Data, has been demultiplexed taken into account that some sequences
may be also in reverse complementary orientation ('--revcomp' setting).
Sequences where reverse complementary indexes have been found 
were reverse complemented, so all the sequences are in uniform orientation in the files.

Sequence orientation in 'demultiplex_out' reflects the indexes orientation: i.e. 
1) if only single-end indexes have been specified, and these indexes are attached to 3'-end of a sequence,
then sequence orientation is 3'-5'.
2) if only single-end indexes have been specified, and these indexes are attached to 5'-end of a sequence,
then sequence orientation is 5'-3'.
3) if paired-end indexes have been specified (both ends of the sequence were supplemented with indexes),
and indexes in the file were specified as 5'_indexes followed by 3'_indexes (fwd_index...rev_index),
then sequence orientation is 5'-3'.\n

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.\n

Core command -> 
cutadapt $indexes_file_in $error_rate $no_indels --revcomp $overlap $minlen

Summary of sequence counts in 'seq_count_summary.txt'

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#cutadapt v4.4 for demultiplexing
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0 for validating indexes file
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

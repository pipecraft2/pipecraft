#!/bin/bash

# Demultiplex PAIRED-END reads
# Demultiplexing of paired-end reads in mixed orientation using single-end or paired-end indexes is supported.
# Input = a directory with fastq/fasta files (R1.fastq; R2.fastq); and indexes file in fasta format (header as a sample name).

##########################################################
###Third-party applications:
#cutadapt v3.5
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #Distributed under the MIT license
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
#pigz v2.4
#python3 with biopython
##################################################################

#Load variables
regex='[^/]*$'
oligos_file_path=$(echo $index_file | grep -oP "$regex")
oligos_file=$(basename $oligos_file_path) #basename, needed for macOS
indexes_file=$(printf "/extraFiles/$oligos_file")
extension=$fileFormat
error_rate="-e ${index_mismatch}"
search_window=${search_window}
if [ "$no_indels" = true ] ; then
    no_indels=$"--no-indels"
else
    no_indels=''
fi
minlen=$"--minimum-length ${min_seq_length}"
cores=$"--cores ${cores}"
overlap=$"--overlap ${overlap}"

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
#python module for assigning sample names as based on indexes file
run_python_module=$"python3 /scripts/submodules/assign_sample_names.demuxModule.py $indexes_file"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env
### Check barcodes file
check_indexes_file

### Process file
printf "Checking files ...\n"
while read LINE; do

    #Write file name without extension
    inputR1=$(echo $LINE | sed -e "s/.$extension//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')

    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        #$extension will be $newextension
    check_gz_zip_PE
    ### Check input formats (fastq/fasta supported)
    check_extension_fastx

    ### Check if dual indexes or single indexes and prepare workflow accordingly
    if [[ $tag == "dual" ]]; then
        #dual indexes
        #make rev barcodes file
        seqkit seq --quiet -n tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
        sed -e 's/^/>/' > tempdir2/sample_names.txt
        grep "\..." tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
        awk 'BEGIN{FS="."}{print $4}' > tempdir2/index_rev.temp
        touch tempdir2/index_rev.fasta
        i=1
        p=$"p"
        while read HEADER; do
            echo $HEADER >> tempdir2/index_rev.fasta
            sed -n $i$p tempdir2/index_rev.temp >> tempdir2/index_rev.fasta
            i=$((i + 1))
        done < tempdir2/sample_names.txt
        rm tempdir2/index_rev.temp
        #make fwd barcodes file
        sed -e 's/\.\.\..*//' < tempdir2/ValidatedBarcodesFileForDemux.fasta.temp > tempdir2/index_fwd.fasta

        #unique index names for assigning sample names on demux files after cutadapt demux
        seqkit rmdup --quiet tempdir2/index_fwd.fasta --by-seq -w 0 > tempdir2/index_fwd.uniq.fasta
        seqkit replace --quiet tempdir2/index_fwd.uniq.fasta -w 0 -p .+ -r "indexF_{nr}" > tempdir2/index_fwd.uniq.renamed.fasta
        seqkit rmdup --quiet tempdir2/index_rev.fasta --by-seq -w 0 > tempdir2/index_rev.uniq.fasta
        seqkit replace --quiet tempdir2/index_rev.uniq.fasta -w 0 -p .+ -r "indexR_{nr}" > tempdir2/index_rev.uniq.renamed.fasta

        #add search window size to indexes
        sed -e '/^>/!s/^/search_window/' tempdir2/index_fwd.uniq.renamed.fasta > $output_dir/index_fwd.fasta
        sed -i "s/search_window/XN{$search_window}/" $output_dir/index_fwd.fasta
        sed -e '/^>/!s/^/search_window/' tempdir2/index_rev.uniq.renamed.fasta> $output_dir/index_rev.fasta
        sed -i "s/search_window/XN{$search_window}/" $output_dir/index_rev.fasta

        #assign demux variables
        indexes_file_round1=$"-g file:$output_dir/index_fwd.fasta -G file:$output_dir/index_rev.fasta"
        indexes_file_round2=$"-g file:index_fwd.fasta -G file:index_rev.fasta"
        outR1=$"-o $output_dir/round1-{name1}-{name2}.R1.$newextension"
        outR2=$"-p $output_dir/round1-{name1}-{name2}.R2.$newextension"
        outR2_round2=$"-o round2-{name1}-{name2}.R2.$newextension"
        outR1_round2=$"-p round2-{name1}-{name2}.R1.$newextension"
        input_for_round2_R1=$"round1-unknown-unknown.R1"
        input_for_round2_R2=$"round1-unknown-unknown.R2"
    else
        #single indexes
        # Add search window size to indexes 
        sed -i '/^>/!s/^/search_window/' tempdir2/ValidatedBarcodesFileForDemux.fasta.temp 
        sed -i "s/search_window/XN{$search_window}/" tempdir2/ValidatedBarcodesFileForDemux.fasta.temp 
        #Move edited indexes file with window size into output_dir 
        mv tempdir2/ValidatedBarcodesFileForDemux.fasta.temp tempdir2/index_file.fasta
        mv tempdir2/index_file.fasta $output_dir

        #assign demux variables
        indexes_file_round1=$"-g file:$output_dir/index_file.fasta"
        indexes_file_round2=$"-g file:index_file.fasta"
        outR1=$"-o $output_dir/round1-{name}.R1.$newextension"
        outR2=$"-p $output_dir/round1-{name}.R2.$newextension"
        outR2_round2=$"-o round2-{name}.R2.$newextension"
        outR1_round2=$"-p round2-{name}.R1.$newextension"
        input_for_round2_R1=$"round1-unknown.R1"
        input_for_round2_R2=$"round1-unknown.R2"
    fi


    ############################
    ### Start demultiplexing ###
    ############################
    printf "\n# Demultiplexing with $tag indexes ... \n"
    printf "   (this may take some time for large files)\n"

    ### Round1 demux
    checkerror=$(cutadapt --quiet \
    $indexes_file_round1 \
    $error_rate \
    $no_indels \
    $overlap \
    $minlen \
    $cores \
    $outR1 \
    $outR2 \
    $inputR1.$newextension $inputR2.$newextension 2>&1)
    check_app_error

    #Round2 demux (RC; R1 and R2 position switched!)
    cd $output_dir
    checkerror=$(cutadapt --quiet \
    $indexes_file_round2 \
    $error_rate \
    $no_indels \
    $overlap \
    $minlen \
    $cores \
    $outR2_round2 \
    $outR1_round2 \
    $input_for_round2_R2.$newextension $input_for_round2_R1.$newextension 2>&1)
    check_app_error
    
    #Remove round1 unknowns and remane final unknowns
    rm $input_for_round2_R2.$newextension
    rm $input_for_round2_R1.$newextension
    if [[ -f round2-unknown.R1.$newextension ]]; then
        mv round2-unknown.R1.$newextension unknown.R1.$newextension
        mv round2-unknown.R2.$newextension unknown.R2.$newextension
    elif [[ -f round2-unknown-unknown.R1.$newextension ]]; then
        mv round2-unknown-unknown.R1.$newextension unknown.R1.$newextension
        mv round2-unknown-unknown.R2.$newextension unknown.R2.$newextension
    fi

    ### Merge demux outputs from round1 and round2
    ls | grep "R1.$newextension" | grep "round1-" > demux_R1_list.txt #ls | grep "R1.$newextension" | grep "round1-" | grep -v "unknown" > demux_R1_list.txt
    while read DEMUXFILES; do
        R1=$(echo $DEMUXFILES | sed -e 's/round1-//' )
        R2=$(echo $DEMUXFILES | sed -e 's/round1-//' | sed -e 's/R1/R2/' )
        if [[ -f round2-$R1 ]]; then
            cat round1-$R1 round2-$R1 > $R1
            rm round1-$R1
            rm round2-$R1
        else
            echo "round2-$R1 does not exist"
        fi
        if [[ -f round2-$R2 ]]; then
            cat round1-$R2 round2-$R2 > $R2
            rm round1-$R2
            rm round2-$R2
        else
            echo "round2-$R2 does not exist"
        fi
    done < demux_R1_list.txt && rm demux_R1_list.txt
    cd ..
done < tempdir2/paired_end_files.txt

# Make patterns file for assigning sample names to files if using DUAL indexes
if [[ $tag = "dual" ]]; then
    ls $output_dir | grep ".R1.$newextension" | \
    sed -e "s/\.R1\.$newextension//" | \
    grep -E -v "round2|round1|unknown" > tempdir2/demux_files.txt

    #Assign sample names
    export newextension
    $run_python_module

    # Move un-named combination fastq files to separate folder
    cd $output_dir
    for unnamed in indexF_*-indexR_*; do
        if [[ -f $unnamed ]]; then
            mkdir -p unnamed_index_combinations
            mv $unnamed unnamed_index_combinations
        fi
    done
    for unnamed in unknown-indexR_*; do
        if [[ -f $unnamed ]]; then
            mkdir -p unnamed_index_combinations
            mv $unnamed unnamed_index_combinations
        fi
    done
    for unnamed in indexF_*-unknown*; do
        if [[ -f $unnamed ]]; then
            mkdir -p unnamed_index_combinations
            mv $unnamed unnamed_index_combinations
        fi
    done
    # Delete empty files in unnamed_index_combinations
    find unnamed_index_combinations -empty -type f -delete
    cd ..

    cp tempdir2/index_fwd.uniq.renamed.fasta $output_dir/unnamed_index_combinations && mv $output_dir/unnamed_index_combinations/index_fwd.uniq.renamed.fasta $output_dir/unnamed_index_combinations/index_fwd.fasta
    cp tempdir2/index_rev.uniq.renamed.fasta $output_dir/unnamed_index_combinations && mv $output_dir/unnamed_index_combinations/index_rev.uniq.renamed.fasta $output_dir/unnamed_index_combinations/index_rev.fasta

fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats_demux

# Add seq count in unnamed_index_combinations to seq_count_summary.txt
cd $output_dir
if [[ -d unnamed_index_combinations ]]; then
  unnamed_index_combinations_count=$(cat unnamed_index_combinations/*.R1.$newextension | seqkit stats --threads 6 -T  | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $4}')
  printf "\n Number of sequences in 'unnamed_index_combinations' dir (*.R1.$newextension): $unnamed_index_combinations_count" >> seq_count_summary.txt
fi
cd ..

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
printf "Files in 'demultiplex_out' directory represent per sample sequence files, 
that were generated based on the specified indexes file ($oligos_file).
[demultiplex_out/index_*fasta files(s) is $oligos_file but with added search window size for cutadapt].

Paired-end data, has been demultiplexed taken into account that some sequences
may be also in reverse complementary orientation.\n
Output R1 and R2 reads are synchronized for merging paired-end data. 

Files in 'unnamed_index_combinations' directory [if present; only when using dual indexes] represent
index combinations that do not correspond to combinations used in indexes file ($oligos_file).

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.
RUNNING THE PROCESS SEVERAL TIMES IN THE SAME DIRECTORY WILL OVERWRITE ALL THE OUTPUTS!

Summary of sequence counts in 'seq_count_summary.txt'

Total run time was $runtime sec.

##########################################################
###Third-party applications used for this process [PLEASE CITE]:
#cutadapt v3.5 for demultiplexing
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.0.0 for validating indexes file and adjusting sample names
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################" > $output_dir/README.txt

###Done, files in $output_dir folder
printf "\nDONE\n"
printf "Data in directory $output_dir\n"
printf "Summary of sequence counts in 'seq_count_summary.txt'\n"
printf "Check README.txt file in $output_dir for further information about the process.\n\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$newextension"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"
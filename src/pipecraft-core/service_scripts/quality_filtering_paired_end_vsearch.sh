#!/bin/bash

# Quality filter PAIRED-END sequencing data with vsearch
# Input = paired-end fastq files

############################
###Third-party applications:
#vsearch v2.23.0
#seqkit v2.3.0
#pigz v2.4
############################

# Checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"

#load variables
export fileFormat 
maxee=$"--fastq_maxee ${maxee}"
maxns=$"--fastq_maxns ${maxNs}"
minlen=$"--fastq_minlen ${min_length}"
cores=$"--threads ${cores}"
qmax=$"--fastq_qmax ${qmax}"
qmin=$"--fastq_qmin ${qmin}"
trunc_length=${trunc_length}
max_length=${max_length}
maxee_rate=${maxee_rate}
truncqual=${truncqual}
truncee=${truncee}

#Source for functions
source /scripts/submodules/framework.functions.sh

#additional options, if selection != undefined
if [[ $max_length == null ]] || [[ -z $max_length ]]; then
    max_length=$""
else
    max_length=$"--fastq_maxlen $max_length"
fi
if [[ $maxee_rate == null ]] || [[ -z $maxee_rate ]]; then
    maxee_rate=$""
else
    maxee_rate=$"--fastq_maxee_rate $maxee_rate"
fi
if [[ $trunc_length == null ]] || [[ -z $trunc_length ]]; then
    trunc_length=$""
else
    trunc_length=$"--fastq_trunclen $trunc_length"
fi
if [[ $truncqual == null ]] || [[ -z $truncqual ]]; then
    truncqual=$""
else
    truncqual=$"--fastq_truncqual $truncqual"
fi
if [[ $truncee == null ]] || [[ -z $truncee ]]; then
    truncee=$""
else
    truncee=$"--fastq_truncee $truncee"
fi

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "vsearch paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = quality filtering"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 3 -mindepth 1 -type d | grep "assembled_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = quality filtering"
    DIRS=$(pwd)
    printf "\n workingDir = $DIRS \n"
fi

#############################
### Start of the workflow ###
#############################
### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)

    cd $seqrun
    if [[ $multiDir == "TRUE" ]]; then
        ### Check if the dir has the specified file extension; if not then ERROR
        count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
        if [[ $count == 0 ]]; then
            printf '%s\n' "ERROR]: cannot find files with specified extension '$fileFormat' in dir $seqrun.
            Please check the extension of your files and specify again.
            >Quitting" >&2
            end_process
        fi
        #output dir
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/qualFiltered_out"
        export output_dir
       
        ### Prepare working env and check paired-end data
        prepare_PE_env
    else
        #output dir
        output_dir=$"/input/qualFiltered_out"
        export output_dir
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check paired-end data
        prepare_PE_env
    fi

    ### Process samples
    while read LINE; do
        #Read in R1 and R2 file names; without extension
        inputR1=$(echo $LINE | sed -e "s/.$fileFormat//")
        inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
        ## Preparing files for the process
        printf "\n____________________________________\n"
        printf "Processing $inputR1 and $inputR2 ...\n"
        #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        check_gz_zip_PE
        ### Check input formats (fastq supported)
        check_extension_fastq

        ###############################
        ### Start quality filtering ###
        ###############################
        mkdir -p tempdir
    
        printf "vsearch --fastq_filter \
        $inputR1.$extension \
        $maxee \
        $maxns \
        $trunc_length \
        $minlen \
        $cores \
        $qmax \
        $qmin \
        $max_length \
        $maxee_rate \
        $truncqual \
        $truncee \
        --fastqout tempdir/$inputR1.$extension"

        #R1
        checkerror=$(vsearch --fastq_filter \
        $inputR1.$extension \
        $maxee \
        $maxns \
        $trunc_length \
        $minlen \
        $cores \
        $qmax \
        $qmin \
        $max_length \
        $maxee_rate \
        $truncqual \
        $truncee \
        --fastqout tempdir/$inputR1.$extension 2>&1)
        check_app_error

        #R2
        checkerror=$(vsearch --fastq_filter \
        $inputR2.$extension \
        $maxee \
        $maxns \
        $trunc_length \
        $minlen \
        $cores \
        $qmax \
        $qmin \
        $max_length \
        $maxee_rate \
        $truncqual \
        $truncee \
        --fastqout tempdir/$inputR2.$extension 2>&1)
        check_app_error

        #If outputs are not empty, then synchronize R1 and R2
        if [ -s tempdir/$inputR1.$extension ]; then
            if [ -s tempdir/$inputR2.$extension ]; then
                printf "\nSynchronizing R1 and R2 reads (matching order for paired-end reads merging)\n"
                cd tempdir
                checkerror=$(seqkit pair -1 $inputR1.$extension -2 $inputR2.$extension 2>&1)
                check_app_error

                rm $inputR1.$extension
                rm $inputR2.$extension
                mv $inputR1.paired.$extension $inputR1.$extension
                mv $inputR2.paired.$extension $inputR2.$extension
                cd ..

                #Move final files to $output_dir
                mv tempdir/$inputR1.$extension $output_dir/$inputR1.$extension
                mv tempdir/$inputR2.$extension $output_dir/$inputR2.$extension

            else
                printf '%s\n' "WARNING]: $inputR2 has 0 seqs after filtering (no output for that sample)"
                rm tempdir/$inputR1.$extension
                rm tempdir/$inputR2.$extension
            fi
        else
            printf '%s\n' "WARNING]: $inputR1 has 0 seqs after filtering (no output for that sample)"
            rm tempdir/$inputR1.$extension
            rm tempdir/$inputR2.$extension
        fi
    done < tempdir2/paired_end_files.txt

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    printf "\nCleaning up and compiling final stats files ...\n"
    clean_and_make_stats
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file
    printf "# Quality filtering was performed using vsearch (see 'Core command' below for the used settings).
    Start time: $start_time
    End time: $(date)
    Runtime: $runtime seconds

    Files in 'qualFiltered_out':
    ----------------------------
    # *.$extension          = quality filtered sequences in FASTQ format.
    # seq_count_summary.txt = summary of sequence counts per sample.

    Core commands -> 
    quality filtering: vsearch --fastq_filter input_file $maxee $maxns $trunc_length $minlen $cores $qmax $qmin $max_length $maxee_rate $truncqual $truncee --fastqout output_file
    Synchronizing R1 and R2 reads (matching pair mates): seqkit pair -1 inputR1 -2 inputR2

    ##############################################
    ###Third-party applications for this process:
    #vsearch (version $vsearch_version) for quality filtering
        #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
        #https://github.com/torognes/vsearch
    #seqkit (version $seqkit_version) for synchronizing R1 and R2 after filtering
        #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
        #https://bioinf.shenwei.me/seqkit/
    ##############################################" > $output_dir/README.txt

    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
else
    echo "workingDir=$output_dir"
fi
echo "fileFormat=$extension"
echo "readType=paired_end"

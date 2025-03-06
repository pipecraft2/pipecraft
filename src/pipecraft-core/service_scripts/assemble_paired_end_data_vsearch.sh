#!/bin/bash

# ASSEMBLE PAIRED-END data with vsearch
# Input = paired-end fastq files. Paired-end data identifiers --> -R1 | _R1 | .R1

################################################
###Third-party applications:
#vsearch, pigz v2.4
################################################

# checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# vsearch (version $vsearch_version)\n"
printf "# pipeline = $pipeline\n"
printf "# service = $service\n"

# load variables
fastq_minoverlen="--fastq_minovlen ${min_overlap}"
fastq_minmergelen="--fastq_minmergelen ${min_length}"
fastq_allowmergestagger=${allow_merge_stagger}
include_R1=$include_only_R1
fastq_maxdiffs="--fastq_maxdiffs ${max_diffs}"
fastq_maxns="--fastq_maxns ${max_Ns}"
fastq_maxmergelen="--fastq_maxmergelen ${max_length}"
fastq_qmax=$fastq_qmax
notmerged_files=$keep_disjointed

# Source for functions
source /scripts/submodules/framework.functions.sh

## check if working with multiple runs or with a single sequencing run
# if working with multiRunDir, and the previous step was CUT PRIMERS
if [[ -f "$workingDir/.prev_step.temp" ]]; then
    prev_step=$(cat $workingDir/.prev_step.temp) # for checking previous step (output from cut_primers_paired_end_reads.sh)
    printf "# prev_step = $prev_step\n"
fi
if [[ -d "/input/multiRunDir" ]] && [[ $pipeline == "vsearch_OTUs" || $pipeline == "UNOISE_ASVs" ]] && [[ $prev_step == "cut_primers" ]]; then
    echo "vsearch paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = assembling paired-end reads (after cut primers)"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "primersCut_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
    rm $workingDir/.prev_step.temp
# if working with multiRunDir, but the previous step was not CUT PRIMERS
elif [[ -d "/input/multiRunDir" ]] && [[ $prev_step != "cut_primers" ]]; then
    echo "vsearch paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = assembling paired-end reads"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*" [replacing ^./ with sed]
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
# if working with individual run (within Quick Tools or pipeline)
else
    echo "Working with individual sequencing run"
    echo "Process = assembling paired-end reads"
    DIRS=$(pwd)
    printf "\n workingDir = $DIRS \n"
fi

### Check file formatting for FASTQ 
if [[ $fileFormat == "fastq" ]] || [[ $fileFormat == "fq" ]] || [[ $fileFormat == "fastq.gz" ]] || [[ $fileFormat == "fq.gz" ]]; then
    :
else
    printf '%s\n' "ERROR]: $fileFormat formatting not supported here!
Supported extensions: fastq(.gz), fq(.gz).
>Quitting" >&2
    end_process
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
        workingDir="/input/multiRunDir/$seqrun"
        export workingDir # real WD for R
        #output dir
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/assembled_out"
        export output_dir
        
        ### Check if the dir has the specified file extension; if not then ERROR
        count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
        if [[ $count == 0 ]]; then
            printf '%s\n' "ERROR]: cannot find files with specified extension '$fileFormat' in dir $seqrun.
            Please check the extension of your files and specify again.
            >Quitting" >&2
            end_process
        fi
    # if working with individual run (within Quick Tools or pipeline)
    else
        workingDir=$seqrun
        export workingDir
        #output dir
        output_dir=$"/input/assembled_out"
        export output_dir
        # Check if files with specified extension exist in the dir
            # check only at the denoising step, not merging step
        if [[ $workingDir == "/input" ]]; then
            first_file_check
            prepare_PE_env
        fi
    fi

    ### Check if files with specified extension exist in the dir
    first_file_check
    ### Prepare working env and check paired-end data
    prepare_PE_env
    ### Process samples
    while read LINE; do
        #Read in R1 and R2 file names; without extension
        inputR1=$(echo $LINE | sed -e "s/.$fileFormat//")
        inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
        ## Preparing files for the process
        printf "\n____________________________________\n"
        printf "Checking $inputR1 and $inputR2 ...\n"
        #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        check_gz_zip_PE
        ### Check input formats (fastq supported)
        check_extension_fastq

        ########################
        ### Start assembling ###
        ########################
        fastqout=$(echo $inputR1 | sed -E "s/$read_R1.*//")

        #variables for not_merged output files
        if [[ $notmerged_files == "TRUE" ]] || [[ $notmerged_files == "true" ]]; then
            mkdir -p $output_dir/not_assembled_paired_end_reads
            fastqout_notmerged_fwd="--fastqout_notmerged_fwd $output_dir/not_assembled_paired_end_reads/$inputR1.notAssembled.$extension"
            fastqout_notmerged_rev="--fastqout_notmerged_rev $output_dir/not_assembled_paired_end_reads/$inputR2.notAssembled.$extension"
        fi
        if [[ $fastq_allowmergestagger == "TRUE" ]] || [[ $fastq_allowmergestagger == "true" ]]; then
            allowmergestagger=$"--fastq_allowmergestagger"
        fi 
        #When including R1 to the assembled output, then include fastqout_notmerged_fwd (in case notmerged_files=FALSE)
        if [[ $include_R1 == "TRUE" ]] || [[ $include_R1 == "true" ]]; then
            if [[ $notmerged_files == "FALSE" ]] || [[ $notmerged_files == "false" ]]; then
                mkdir -p $output_dir/not_assembled_paired_end_reads
                fastqout_notmerged_fwd="--fastqout_notmerged_fwd $output_dir/not_assembled_paired_end_reads/$inputR1.notAssembled.$extension"
            fi
        fi

        checkerror=$(vsearch --quiet --fastq_mergepairs $inputR1.$extension \
        --reverse $inputR2.$extension \
        $fastq_minoverlen \
        $fastq_minmergelen \
        $allowmergestagger \
        $fastq_maxdiffs \
        $fastq_maxns \
        $fastq_maxmergelen \
        $fastqout_notmerged_fwd \
        $fastqout_notmerged_rev \
        --fastq_qmax $fastq_qmax \
        --fastq_qmaxout $fastq_qmax \
        --fastqout $output_dir/$fastqout.$extension 2>&1)
        check_app_error

        #Include R1 reads to assembled data set if include_R1 = TRUE
        if [[ $include_R1 == "TRUE" ]] || [[ $include_R1 == "true" ]]; then
            printf '%s\n' "include only R1 = TRUE, including also unmerged R1 reads to the assembled output"
            cat $output_dir/not_assembled_paired_end_reads/$inputR1.notAssembled.$extension >> $output_dir/$fastqout.$extension
        fi
        ### Check if assembled output is empty; if yes, then report WARNING
        if [ -s $output_dir/$fastqout.$extension ]; then
            :
        else
            printf '%s\n' "WARNING]: after assembling, $fastqout.$extension has 0 seqs (no output)" >&2
        fi
    done < tempdir2/paired_end_files.txt

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    printf "\nCleaning up and compiling final stats files ...\n"
    clean_and_make_stats_assemble
    if [[ $include_R1 == "TRUE" ]] || [[ $include_R1 == "true" ]]; then
        if [[ $notmerged_files == "FALSE" ]] || [[ $notmerged_files == "false" ]]; then
            rm -r $output_dir/not_assembled_paired_end_reads
        fi
    fi
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file for demultiplexed reads
    printf "# Paired-end data was assembled using vsearch (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'assembled_out' directory represent assembled paired-end files.

If include only R1 = TRUE, then the unassembled R1 reads have been added to the set of assembled reads per sample.
This may be relevant when working with e.g. ITS2 sequences, because ITS2 region in some taxa is too long for assembly, 
therefore discarded completely after assembly process. Thus, including also unassembled R1 reads, partial ITS2 sequences 
for these taxa will be represented in the final output. 
If include only R1 option = TRUE, then other specified options (length, max error rate etc.) have not been 
applied to R1 reads in the 'assembled' file. Thus, additional quality filtering (if this was done before assembling) 
should be run on the 'assembled' data.\n
NOTE RUNNING THE PROCESS SEVERAL TIMES IN THE SAME DIRECTORY WILL OVERWRITE ALL THE OUTPUTS!

Core command -> 
vsearch --fastq_mergepairs input.R1 --reverse input.R2 $fastq_minoverlen $fastq_minmergelen $allowmergestagger $fastq_maxdiffs $fastq_maxns $fastq_maxmergelen $fastqout_notmerged_fwd $fastqout_notmerged_rev --fastq_qmax $fastq_qmax --fastq_qmaxout $fastq_qmax --fastqout output_file

Summary of sequence counts in 'seq_count_summary.txt'

############################################
###Third-party applications for this process:
#vsearch (version $vsearch_version) for assembling paired-end reads
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
############################################" > $output_dir/README.txt

    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done

###Done, files in $output_dir folder
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
    # var for multiRunDir pipe
    printf "merge_reads" > $workingDir/.prev_step.temp
else
    echo "workingDir=$output_dir"
fi
echo "fileFormat=$extension"
echo "readType=single_end"

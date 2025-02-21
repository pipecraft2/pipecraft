#!/bin/bash

# Quality filter PAIRED-END sequencing data with dada2
# Input = paired-end fastq files

############################
###Third-party applications:
# dada2 v1.28
# seqkit v2.3.0
############################

# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# DADA2 version: $dada2_version\n"
printf "# seqkit version: $seqkit_version\n"

#load env variables
readType=${readType}
dataFormat=${dataFormat}
workingDir=${workingDir}

### variables
# maxEE       = discard sequences with more than the specified number of expected errors
# maxN        = discard sequences with more than the specified number of Ns (ambiguous bases)
# truncQ      = truncate reads at the first instance of a quality score less than or equal to truncQ
# truncLen    = truncate reads after truncLen bases (applies to R1 reads when working with paired-end data)
# truncLen_R2 = truncate R2 reads after truncLen bases
# minLen      = remove reads with length less than minLen. minLen is enforced after all other trimming and truncation
# maxLen      = remove reads with length greater than maxLen. maxLen is enforced on the raw reads
# minQ        = after truncation, reads contain a quality score below minQ will be discarded
# matchIDs    = If true, then double-checking (with seqkit pair) that only paired reads that share ids are outputted
                    #WORK WITH SEQKIT for matchIDs = true, because sometimes DADA2 CANNOT automatically identify paired-end headers

#Source for functions
source /scripts/submodules/framework.functions.sh

## check if working with multiple runs or with a single sequencing run
# if working with multiRunDir, and the previous step was CUT PRIMERS
if [[ -f "$workingDir/.prev_step.temp" ]]; then
    prev_step=$(cat $workingDir/.prev_step.temp) # for checking previous step (output from cut_primers_paired_end_reads.sh)
fi
if [[ $pipeline == "DADA2_ASVs" ]] && [[ $prev_step == "cut_primers" ]]; then
    echo "DADA2 paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = quality filtering (after cut primers)"
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
    echo "DADA2 paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = quality filtering"
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
    echo "Process = quality filtering"
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

    ### Process samples with dada2 filterAndTrim function in R
    printf "# Running DADA2 filterAndTrim in $seqrun \n"
    Rlog=$(Rscript /scripts/submodules/dada2_PE_filterAndTrim.R 2>&1)
    echo $Rlog > $output_dir/dada2_PE_filterAndTrim.log 
    wait
    #format R-log file
    sed -i "s/;; /\n/g" $output_dir/dada2_PE_filterAndTrim.log 


    ### Synchronizing R1 and R2 reads if $matchIDs == "true" - WORK WITH SEQKIT for matchIDs = true, because sometimes DADA2 CANNOT automatically identify paired-end headers
    if [[ $matchIDs == "true" ]] || [[ $matchIDs == "TRUE" ]]; then
        while read LINE; do
            #Read in R1 and R2 file names; without extension
            samp_name=$(basename $LINE | awk -F\\${read_R1} '{print$1}')
            #If outputs are not empty, then synchronize R1 and R2
            if [[ -s $output_dir/$samp_name\_R1.$fileFormat ]]; then
                if [[ -s $output_dir/$samp_name\_R2.$fileFormat ]]; then
                    printf "\nSynchronizing $samp_name R1 and R2 reads\n"
                    cd $output_dir
                    checkerror=$(seqkit pair -1 $samp_name\_R1.$fileFormat -2 $samp_name\_R2.$fileFormat 2>&1)
                    check_app_error

                    rm $samp_name\_R1.$fileFormat
                    rm $samp_name\_R2.$fileFormat
                    mv $samp_name\_R1.paired.$fileFormat $samp_name\_R1.$fileFormat
                    mv $samp_name\_R2.paired.$fileFormat $samp_name\_R2.$fileFormat
                    cd ..
                fi
            else
                printf "NOTE: all reads descarded from $samp_name\n"
            fi
        done < tempdir2/paired_end_files.txt
    fi

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    if [[ $debugger != "true" ]]; then
        if [[ -d "tempdir2" ]]; then
            rm -rf tempdir2
        fi
        if [[ -d "primersCut_out/tempdir2" ]]; then
            rm -rf primersCut_out/tempdir2
        fi
        #rm $output_dir/dada2_PE_filterAndTrim.log
    fi
    ### end pipe if no outputs were generated
    outfile_check=$(ls $output_dir/*.$fileFormat 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then 
        :
    else 
        printf '%s\n' "ERROR]: no output files generated after quality filtering ($output_dir). Adjust settings or check sample identifier 'read_R1/R2' so that all sample names would be unique.
        >Quitting" >&2
        end_process
    fi

    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file
    printf "# Quality filtering was performed using dada2 (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'qualFiltered_out':
----------------------------
    # *.$fileFormat            = quality filtered sequences per sample.
    # seq_count_summary.csv    = summary of sequence counts per sample.
    # *.rds                    = R objects for dada2.

Core command -> 
filterAndTrim(inputR1, outputR1, inputR2, outputR2, maxN = $maxN, maxEE = c($maxEE, $maxEE), truncQ = $truncQ, truncLen = c($truncLen, $truncLen_R2), maxLen = $maxLen, minLen = $minLen, minQ=$minQ, rm.phix = TRUE)

Total run time was $runtime sec.
##############################################
###Third-party applications for this process:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
# seqkit (version $seqkit_version) for synchronizing R1 and R2 after filtering (when matchIDs = TRUE)
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

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
else
    echo "workingDir=$output_dir"
fi
echo "fileFormat=$fileFormat"
echo "readType=paired_end"


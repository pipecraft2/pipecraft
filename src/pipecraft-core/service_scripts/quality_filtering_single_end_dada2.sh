#!/bin/bash

# Quality filter SINGLE-END sequencing data with dada2
# Input = single-end fastq files

################################################
###Third-party applications:
#dada2 v1.28
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# DADA2 version: $dada2_version\n"

#load env variables
readType=${readType}
extension=${fileFormat}
workingDir=${workingDir}

#Source for functions
source /scripts/submodules/framework.functions.sh

## check if working with multiple runs or with a single sequencing run
# if working with multiRunDir, and the previous step was CUT PRIMERS
if [[ -f "$workingDir/prev_step.temp" ]]; then
    prev_step=$(cat $workingDir/prev_step.temp) # for checking previous step (output from cut_primers_paired_end_reads.sh)
fi
if [[ $pipeline == "DADA2_ASVs" ]] && [[ $prev_step == "cut_SE_primers" ]]; then
    echo "DADA2 single-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = quality filtering (after cut primers)"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "primersCut_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
    rm $workingDir/prev_step.temp
# if working with multiRunDir, but the previous step was not CUT PRIMERS
elif [[ -d "/input/multiRunDir" ]] && [[ $prev_step != "cut_SE_primers" ]]; then
    echo "DADA2 single-end pipeline with multiple sequencing runs in multiRunDir"
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
        ### Prepare working env and check single-end data
        prepare_SE_env
    else
        #output dir
        output_dir=$"/input/qualFiltered_out"
        export output_dir
        # Check if files with specified extension exist in the dir
        first_file_check
        ### Prepare working env and check single-end data
        prepare_SE_env
    fi

    ### Process samples with dada2 filterAndTrim function in R
    printf "# Running DADA2 filterAndTrim in $seqrun\n"
    Rlog=$(Rscript /scripts/submodules/dada2_SE_filterAndTrim.R 2>&1)
    echo $Rlog > $output_dir/qFilt.log 
    wait
    printf "\n DADA2 filterAndTrim completed \n"

    #format R-log file
    sed -i "s/;; /\n/g" $output_dir/qFilt.log 

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        #rm $output_dir/qFilt.log 
    fi

    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file
    printf "# Quality filtering was performed using dada2 (see 'Core command' below for the used settings).

Files in 'qualFiltered_out':
----------------------------
# *.$extension             = quality filtered sequences per sample.
# seq_count_summary.csv    = summary of sequence counts per sample.
# *.rds                    = R objects for dada2.

Core command -> 
filterAndTrim(inputR1, outputR1, maxN = $maxN, maxEE = $maxEE, truncQ = $truncQ, truncLen = $truncLen, maxLen = $maxLen, minLen = $minLen, minQ = $minQ, rm.phix = TRUE)

Total run time was $runtime sec.
##############################################
###Third-party applications for this process:
#dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
##############################################" > $output_dir/README.txt

    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
else
    echo "workingDir=$output_dir"
fi
echo "fileFormat=$extension"
echo "readType=single_end"


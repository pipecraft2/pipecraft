#!/bin/bash

# Chimera filtering with DADA2 removeBimeraDenovo function for ASVs workflow [MIXED orient amplicons].

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
workingDir=${workingDir}

## settings
# method = list ["consensus", "pooled", "per-sample"]. 

#Source for functions
source /scripts/submodules/framework.functions.sh

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "DADA2 paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = chimera filtering"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")
    echo "Working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    cd /input
    echo "Working with individual sequencing run"
    echo "Process = chimera filtering (for MIXED amplicons)"
    DIRS=$"/input" # looking those DIRS to merge data: primersCut_out/fwd_orient/denoised_assembled.dada2    /input/primersCut_out/rev_orient/denoised_assembled.dada2"
    printf "\n workingDirs = $DIRS \n"
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
        output_dir1=$"/input/multiRunDir/${seqrun%%/*}/chimeraFiltered_out.dada2"
        output_dir2=$"/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2"
        export output_dir1
        export output_dir2
        #WDs for fwd_orient and rev_orient
        workingDir_fwd=$"/input/multiRunDir/${seqrun%%/*}/primersCut_out/fwd_orient/denoised_assembled.dada2"
        export workingDir_fwd
        workingDir_rev=$"/input/multiRunDir/${seqrun%%/*}/primersCut_out/rev_orient/denoised_assembled.dada2"
        export workingDir_rev

        # check for ASVs_table.denoised.rds
        if [[ ! -f "$workingDir_fwd/ASVs_table.denoised.rds" ]]; then
            printf '%s\n' "ERROR]: cannot find ASVs_table.denoised.rds in $seqrun primersCut_out/fwd_orient/denoised_assembled.dada2.
            Denoising failed? Did not produce any output. Please report this ERROR.
            >Quitting" >&2
            end_process
        fi
        if [[ ! -f "$workingDir_rev/ASVs_table.denoised.rds" ]]; then
            printf '%s\n' "ERROR]: cannot find ASVs_table.denoised.rds in $seqrun primersCut_out/rev_orient/denoised_assembled.dada2.
            Denoising failed? Did not produce any output. Please report this ERROR.
            >Quitting" >&2
            end_process
        fi
    else
        workingDir=$"/input"
        export workingDir
        #output dir
        output_dir1=$"/input/chimeraFiltered_out.dada2"
        output_dir2=$"/input/ASVs_out.dada2"
        export output_dir1
        export output_dir2
        #WDs for fwd_orient and rev_orient
        workingDir_fwd="/input/primersCut_out/fwd_orient/denoised_assembled.dada2"
        export workingDir_fwd
        workingDir_rev="/input/primersCut_out/rev_orient/denoised_assembled.dada2"
        export workingDir_rev

        # check for ASVs_table.denoised.rds
        if [[ ! -f "$workingDir_fwd/ASVs_table.denoised.rds" ]]; then
            printf '%s\n' "ERROR]: cannot find $workingDir_fwd/ASVs_table.denoised.rds.
            Denoising failed? Did not produce any output. Please report this ERROR.
            >Quitting" >&2
            end_process
        fi
        if [[ ! -f "$workingDir_rev/ASVs_table.denoised.rds" ]]; then
            printf '%s\n' "ERROR]: cannot find $workingDir_rev/ASVs_table.denoised.rds.
            Denoising failed? Did not produce any output. Please report this ERROR.
            >Quitting" >&2
            end_process
        fi
    fi
   
    ### Process samples with dada2 removeBimeraDenovo function in R
    printf "# Running DADA2 removeBimeraDenovo for $seqrun\n "
    Rlog=$(Rscript /scripts/submodules/dada2_chimeraFilt_mixed_wf.R 2>&1)
    echo $Rlog > $output_dir1/dada2_chimeraFilt.log
    wait
    #format R-log file
    sed -i "s/;; /\n/g" $output_dir1/dada2_chimeraFilt.log

    # Add "ASV" as a 1st col name
    if [[ -s $output_dir2/ASVs_table.txt ]]; then
        sed -i "1 s|^|ASV|" $output_dir2/ASVs_table.txt
    fi

    ### NOTE: dada2mode = MIXED -> does not paste out chimeric ASVs per-sample

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        #rm $output_dir1/dada2_chimeraFilt.log 
    fi
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file (chimeraFiltered_out.dada2)
    printf "# Chimeras were filtered out with DADA2 removeBimeraDenovo function (method = $method).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'chimeraFiltered_out.dada2':
-----------------------------------
# *.chimFilt_ASVs.fasta = chimera filtered ASVs per sample. 'Size' denotes the abundance of the ASV sequence  
# seq_count_summary.csv = summary of sequence and ASV counts per sample

##############################################
###Third-party applications for this process:
#dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
##############################################" > $output_dir1/README.txt

    #Make README.txt file (ASVs_out.dada2)
    # count features and sequences; outputs variables feature_count, nSeqs, nSample
    count_features "$output_dir2/ASVs_table.txt"

    printf "# ASV table was constructed with DADA2 makeSequenceTable function.

Files in 'ASVs_out.dada2' directory:
-----------------------------------
# ASVs_table.txt                  = ASV-by-sample table (tab delimited file). [denoised and chimera filtered]
# ASVs.fasta                      = FASTA formated representative ASV sequences [denoised and chimera filtered]
# ASVs_table.denoised.nochim.rds  = rds formatted denoised and chimera filtered ASV table (for DADA2)
# ASVs_table.denoised.rds         = rds formatted denoised ASV table (for DADA2)

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample

##############################################
###Third-party applications for this process:
#dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
##############################################" > $output_dir2/README.txt
    
    
    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    else
        cd /input
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
    echo "workingDir=$output_dir2"
fi
echo "fileFormat=fasta"
echo "readType=single_end"



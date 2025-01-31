#!/bin/bash

# Chimera filtering with DADA2 removeBimeraDenovo function for ASVs workflow.

################################################
###Third-party applications:
# dada2 v1.28
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# DADA2 version: $dada2_version\n"

#load env variables
readType=${readType}
dataFormat=${dataFormat}
workingDir=${workingDir}

# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
printf "# DADA2 version: $dada2_version\n"

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
    DIRS=$(find . -maxdepth 3 -mindepth 1 -type d | grep "denoised_assembled.dada2" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"true"
    export multiDir
else
  echo "Working with individual sequencing run"
  cd /input
  DIRS=$(find . -maxdepth 3 -mindepth 1 -type d | grep "denoised_assembled.dada2" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
fi

#############################
### Start of the workflow ###
#############################
### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
declare -a output_feature_tables
declare -a output_fastas

for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)
    cd $seqrun
    
    if [[ $multiDir == "true" ]]; then
        workingDir="/input/multiRunDir/$seqrun"
        export workingDir # real WD for R
        #output dir
        output_dir1=$"/input/multiRunDir/${seqrun%%/*}/chimeraFiltered_out.dada2"
        output_dir2=$"/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2"
        export output_dir1
        export output_dir2
        # Store output files in arrays for multiple runs
        output_feature_tables+=("$output_dir2/ASVs_table.txt")
        output_fastas+=("$output_dir2/ASVs.fasta")
    else
        #output dir
        output_dir1=$"/input/chimeraFiltered_out.dada2"
        output_dir2=$"/input/ASVs_out.dada2"
        export output_dir1
        export output_dir2
        # Store single output files
        output_feature_table="$output_dir2/ASVs_table.txt"
        output_fasta="$output_dir2/ASVs.fasta"
    fi

    # # FOR TESTING: Skip this process if output directory already exists
    if [[ -d $output_dir1 ]]; then
        printf "# skipping chimera filtering\n"
        if [[ $multiDir == "true" ]]; then
            workingDir=$"/input/multiRunDir"
            echo "workingDir=$workingDir"
        else
            echo "workingDir=$output_dir2"
        fi
        echo "fileFormat=fasta"
        echo "readType=single_end"
        exit 0
    fi
   
    ### Process samples with dada2 removeBimeraDenovo function in R
    printf "# Running DADA2 removeBimeraDenovo in $seqrun\n"
    Rlog=$(Rscript /scripts/submodules/dada2_chimeraFilt_wf.R 2>&1)
    echo $Rlog > $output_dir1/dada2_chimeraFilt.log 
    wait
    #format R-log file
    sed -i "s/;; /\n/g" $output_dir1/dada2_chimeraFilt.log 

    # Add "ASV" as a 1st col name
    if [[ -s $output_dir2/ASVs_table.txt ]]; then
        sed -i "1 s|^|ASV|" $output_dir2/ASVs_table.txt
    fi

    ### Compare 'chimera filtered fasta files per sample' and 'NOT chimera 
     # filtered fasta files per sample' to paste out only chimeric sequences per sample
    echo "pasting chimeric seqs"
    path_denoised=$"$workingDir"

    #make dir for chimeras.fasta
    mkdir $output_dir1/chimeras
    #make seqs_count_summary.txt
    touch $output_dir1/chimeras/seq_count_summary.txt
    printf "File\tReads\n" > $output_dir1/chimeras/seq_count_summary.txt

    for chim_filt_file in $output_dir1/*.chimFilt_ASVs.fasta; do
        samp_name=$(basename $chim_filt_file | awk 'BEGIN{FS=".chimFilt_ASVs.fasta"} {print $1}')
        corresponding_denoised_file=$(ls $path_denoised | grep "$samp_name")

        #seqkit paste chimeras
        seqkit grep -w 0 -svf \
        <(seqkit seq -s -w 0 $chim_filt_file) \
        <(seqkit seq -w 0 $path_denoised/$corresponding_denoised_file) \
        > $output_dir1/chimeras/$samp_name.chimeras.fasta

        #delete if chimeras file is empty
        if [[ ! -s $output_dir1/chimeras/$samp_name.chimeras.fasta ]]; then
            printf "$samp_name.chimeras.fasta\t0\n" >> $output_dir1/chimeras/seq_count_summary.txt
            rm $output_dir1/chimeras/$samp_name.chimeras.fasta
        else
            #count chimeric seqs
            seq_count=$(grep -c "^>" $output_dir1/chimeras/$samp_name.chimeras.fasta)
            printf "$samp_name.chimeras.fasta\t$seq_count\n" >> $output_dir1/chimeras/seq_count_summary.txt
        fi
    done

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

Files in 'chimeraFiltered_out.dada2':
-----------------------------------
# *.chimFilt_ASVs.fasta   = chimera filtered ASVs per sample. 'Size' denotes the abundance of the ASV sequence  
# seq_count_summary.txt   = summary of sequence and ASV counts per sample

Files in 'chimeraFiltered_out.dada2/chimeras':
----------------------------------------------
# *.chimeras.fasta        = ASVs per sample that were flagged as chimeras (and thus discarded).

Total run time was $runtime sec.

#################################################
###Third-party applications for this process:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
#################################################" > $output_dir1/README.txt

    #Make README.txt file (ASVs_out.dada2)
    # count features and sequences; outputs variables feature_count, nSeqs, nSample
    count_features "$output_dir2/ASVs_table.txt"

    printf "# ASV table was constructed with DADA2 makeSequenceTable function.

Files in 'ASVs_out.dada2' directory:
------------------------------------
# ASVs_table.txt                  = ASV-by-sample table (tab delimited file). 
# ASVs.fasta                      = FASTA formated representative ASV sequences 
# ASVs_table.denoised.nochim.rds  = rds formatted denoised and chimera filtered ASV table (for DADA2)
# ASVs_table.denoised.nochim.rds  = rds formatted denoised and chimera filtered ASV table (for DADA2)
# ASVs_table.denoised.rds         = rds formatted denoised ASV table (for DADA2) [same .rds is in denoised_assembled.dada2 dir]

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample

##############################################
###Third-party applications for this process:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
##############################################" > $output_dir2/README.txt


    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "true" ]]; then 
        cd /input/multiRunDir
    else
        cd ..
    fi
done

# Validate and combine output files
if [[ $multiDir == "true" ]]; then
    # Check each file in the arrays
    valid_feature_tables=()
    valid_fastas=()
    
    for table in "${output_feature_tables[@]}"; do
        if [[ -f "$table" && -s "$table" ]]; then
            valid_feature_tables+=("$table")
        else
            printf "Warning: Feature table not found or empty: $table\n" >&2
        fi
    done
    
    for fasta in "${output_fastas[@]}"; do
        if [[ -f "$fasta" && -s "$fasta" ]]; then
            valid_fastas+=("$fasta")
        else
            printf "Warning: FASTA file not found or empty: $fasta\n" >&2
        fi
    done
      
    output_feature_table=$(IFS=,; echo "${valid_feature_tables[*]}")
    output_fasta=$(IFS=,; echo "${valid_fastas[*]}")
else
    # Check single run output files
    if [[ ! -f "$output_feature_table" || ! -s "$output_feature_table" ]]; then
        printf "Error: Feature table not found or empty: $output_feature_table\n" >&2
        exit 1
    fi
    if [[ ! -f "$output_fasta" || ! -s "$output_fasta" ]]; then
        printf "Error: FASTA file not found or empty: $output_fasta\n" >&2
        exit 1
    fi
fi

# Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
if [[ $multiDir == "true" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
else
    echo "workingDir=$output_dir2"
fi
echo "fileFormat=fasta"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"




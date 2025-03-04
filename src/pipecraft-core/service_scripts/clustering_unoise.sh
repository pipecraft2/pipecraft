#!/bin/bash

#Input = single-end fasta/fastq files.
#Output = FASTA formated zOTU sequences and zOTU_table.txt, and optionally OTU sequences and OTU_table.txt

# Sequence denoising and clustering

# disabling zOTUs clustering to OTUs for v1.1.0. Number of sequences in OTU table does not match with the seqs in zOTU table (but they should match)


################################################
###Third-party applications:
#vsearch v2.29.4
#GNU Parallel 20210422
#pigz
################################################

# checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"
printf "# pipeline = $pipeline\n"
printf "# service = $service\n"
printf "# clustering dirs = $workingDir\n"

###############################
###############################
#load variables
id=$"--id 1"                                    # positive float (0-1)  # OTU clustering if id < 1. fixed to 1 for v1.1.0 to disable OTU clustering
id_float="1"                                    # fixed to 1 for v1.1.0 to disable OTU clustering
strands=$"--strand ${strands}"                  # both/plus
minsize=$"--minsize ${minsize}"                 # positive integer (default, 8)

#additional options
unoise_alpha=$"--unoise_alpha ${unoise_alpha}"  # positive integer (default, 2)
denoise_level=${denoise_level}                  # list: "global" or "individual"
chimerarm=${remove_chimeras}                    # TRUE or undefined
cores=$"--threads ${cores}"                     # positive integer
abskew=$"--abskew ${abskew}"                    # positive integer (default, 16)
simtype=$"--iddef ${similarity_type}"           # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
maxaccepts=$"--maxaccepts ${maxaccepts}"        # positive integer (default, 1)
maxrejects=$"--maxrejects ${maxrejects}"        # positive integer (default, 32)
mask=$"--qmask ${mask}"                         # list: --qmask dust, --qmask none
###############################
###############################
# Source for functions
source /scripts/submodules/framework.functions.sh

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "UNOISE pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = unoise3"
    cd /input/multiRunDir

    # rm old clustering params if they exist
    if [[ -f "/input/multiRunDir/.clustering_params" ]]; then
        rm /input/multiRunDir/.clustering_params
    fi
    if [[ -f "/input/multiRunDir/.derep_seqs_dirs" ]]; then
        rm /input/multiRunDir/.derep_seqs_dirs
    fi

    # check if prev_step.temp exists
    if [[ -f "$workingDir/.prev_step.temp" ]]; then
        prev_step=$(cat $workingDir/.prev_step.temp) # for checking previous step (output temp file from ITS_extractor.sh)
        printf "# prev_step = $prev_step\n"
        ITS_region_for_clustering=$(cat $workingDir/.ITS_region_for_clustering.temp)
        ITS_full_and_partial=$(cat $workingDir/.ITS_full_and_partial.temp)
        printf "# ITS_region_for_clustering = $ITS_region_for_clustering\n"
        printf "# ITS_full_and_partial = $ITS_full_and_partial (=FALSE if empty)\n"
        rm $workingDir/.prev_step.temp
        rm $workingDir/.ITS_region_for_clustering.temp
        rm $workingDir/.ITS_full_and_partial.temp
    fi
    
    # if working with multiRunDir, and the previous step was ITSx and full_and_partial = FALSE
    if [[ "$prev_step" = "ITSx" ]] && [[ "$ITS_full_and_partial" == "" ]]; then 
        DIRS=$(find . -maxdepth 3 -mindepth 1 -type d | grep "ITSx_out/${ITS_region_for_clustering}" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    # if working with multiRunDir, and the previous step was ITSx and full_and_partial = TRUE
    elif [[ "$prev_step" = "ITSx" ]] && [[ "$ITS_full_and_partial" == "full_and_partial" ]]; then   
        DIRS=$(find . -maxdepth 4 -mindepth 1 -type d | grep "ITSx_out/${ITS_region_for_clustering}/${ITS_full_and_partial}" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    # if working with multiRunDir, and the previous step was not ITSx
    else
        DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "qualFiltered_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    fi
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = unoise3"
    DIRS=$(pwd)
    printf "\n workingDir = $DIRS \n"
fi

#############################
### Start of the workflow ###
#############################
# Store output files in arrays for multiple runs
declare -a output_feature_tables
declare -a output_fastas
 # output files in arrays for if zOTUs are clustered
if [[ $id_float != 1 ]]; then
    declare -a output_feature_tables2
    declare -a output_fastas2
fi

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
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/clustering_out"
        export output_dir
        # Store output files in arrays for multiple runs
        output_feature_tables+=("$output_dir/zOTU_table.txt")
        output_fastas+=("$output_dir/zOTUs.fasta")

        ### Prepare working env and check single-end data
        first_file_check
        prepare_SE_env
    else
        #output dir
        output_dir=$"/input/clustering_out"
        export output_dir
        # Store single output files
        output_feature_table="$output_dir/zOTU_table.txt"
        output_fasta="$output_dir/zOTUs.fasta"
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check single-end data
        prepare_SE_env
    fi

    ### Pre-process samples
    printf "Checking files ... \n"
    for file in *.$fileFormat; do
        #Read file name; without extension
        input=$(echo $file | sed -e "s/.$fileFormat//")
        #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        check_gz_zip_SE
        ### Check input formats
        check_extension_fastx
    done

    #If input is FASTQ then convert to FASTA
    if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
        for file in *.$extension; do
            samp_name=$(basename $file | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
            checkerror=$(seqkit fq2fa -t dna --line-width 0 $file -o $samp_name.fasta 2>&1)
            check_app_error
        done

        was_fastq=$"true"
        export was_fastq
        extension=$"fasta"
        export extension
    fi

    #tempdir
    if [[ -d tempdir ]]; then
        rm -rf tempdir
    fi
    mkdir -p tempdir

    if [[ ! -d dereplicated_sequences ]]; then
        rm -rf dereplicated_sequences
    fi
    mkdir -p dereplicated_sequences

    printf "Dereplication of individual samples ... \n"
    ## Dereplication of individual samples, add sample ID to the header
    derep_rename () {
      samp_name=$(basename $1 | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
      vsearch \
        --derep_fulllength "$1" \
        --relabel_sha1 \
        --output - \
        --fasta_width 0 \
        --sizein --sizeout \
      | sed 's/>.*/&;sample='"$samp_name"'/' > dereplicated_sequences/"$samp_name".fasta
    }
    export -f derep_rename
    find . -maxdepth 1 -name "*.$extension" | parallel -j 1 "derep_rename {}"

    ### Global dereplication
    printf "Dereplicating globally ... \n"
    find dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
    | vsearch \
    --derep_fulllength - \
    --uc tempdir/Glob_derep.uc \
    --output - \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout > $output_dir/Glob_derep.fasta

    ## Denoizing sequences globally
    if [[ $denoise_level == "global" ]]; then
      printf "Denoizing sequences globally ... \n"
    
      ### UNOISE3
      printf "Unoise3 ... \n"
      checkerror=$(vsearch \
      --cluster_unoise $output_dir/Glob_derep.fasta \
      $strands \
      $minsize \
      $unoise_alpha \
      $simtype \
      $mask \
      $maxaccepts \
      $maxrejects \
      $cores \
      --centroids $output_dir/zOTUs.fasta \
      --uc $output_dir/zOTUs.uc \
      --fasta_width 0 \
      --sizein --sizeout 2>&1)
      check_app_error
      
      ## Remove chimera
      printf "Remove chimeras ... \n"

      if [[ $chimerarm == "true" ]]; then
        checkerror=$(vsearch \
        --sortbysize $output_dir/zOTUs.fasta \
        --output - \
        | vsearch --uchime3_denovo - \
        $abskew \
        --nonchimeras $output_dir/zOTUs_noChim.temp.fasta \
        --chimeras $output_dir/UNOISE_Chimeras.fasta 2>&1)
        check_app_error
      
        ## Count number of chimeric sequences
        chimeras=$(grep -c "^>" $output_dir/UNOISE_Chimeras.fasta)

        ## Replace zOTUs with chimera-filtered zOTUs
        rm $output_dir/zOTUs.fasta
        checkerror=$(vsearch --fastx_filter $output_dir/zOTUs_noChim.temp.fasta \
        --fasta_width 0 --fastaout $output_dir/zOTUs.fasta 2>&1)
        check_app_error
        rm $output_dir/zOTUs_noChim.temp.fasta
      fi
    fi  # end of global denoising


    ## Denoizing sequences individually for each sample
    if [[ $denoise_level == "individual" ]]; then
      mkdir -p tempdir_denoize
      mkdir -p tempdir_chimera

      ## Function to denoise and remove chimera for each sample individually 
      denoise_and_chim () {
        
        samp_name=$(basename $1)

        ## Denoise sample
        checkerror=$(vsearch \
        --cluster_unoise "$1" \
          $strands \
          $minsize \
          $unoise_alpha \
          $simtype \
          $mask \
          $maxaccepts \
          $maxrejects \
          --threads 1 \
          --centroids tempdir_denoize/"$samp_name" \
          --fasta_width 0 \
          --sizein --sizeout 2>&1)
        check_app_error
        
        ## Remove chimera
        if [[ $chimerarm == "true" ]]; then
          checkerror=$(vsearch \
            --sortbysize tempdir_denoize/"$samp_name" \
            --output - \
            | vsearch \
            --uchime3_denovo - \
            $abskew \
            --nonchimeras tempdir_chimera/NonChim_"$samp_name" \
            --chimeras tempdir_chimera/Chim_"$samp_name" 2>&1)
          check_app_error
        fi
      }

      export -f denoise_and_chim
      export -f check_app_error
      export -f end_process

      export chimerarm="$chimerarm"
      export strands="$strands"
      export minsize="$minsize"
      export unoise_alpha="$unoise_alpha"
      export simtype="$simtype"
      export mask="$mask"
      export maxaccepts="$maxaccepts"
      export maxrejects="$maxrejects"

      ## Take dereplicated samples and apply denoising function
      printf "Denoizing sequences individually ... \n"
      find dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "denoise_and_chim {}"

      if [[ $chimerarm == "true" ]]; then
        printf "Removing chimeras ... \n"
        find tempdir_chimera -maxdepth 1 -name "Chim_*.fasta" | parallel -j 1 "cat {} >> tempdir/All_chimera.fasta"
        ## Count chimeric sequences
        chimeras=$(grep -c "^>" tempdir/All_chimera.fasta)

        ## Combine and dereplicate denoised sequences
        find tempdir_chimera -maxdepth 1 -name "NonChim_*.fasta" | parallel -j 1 "cat {}" \
        | vsearch \
        --derep_fulllength - \
        --output $output_dir/zOTUs.fasta \
        --uc $output_dir/zOTUs.uc \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout
      
      else
        ## Combine and dereplicate denoised sequences (without chimera removal step)
        find tempdir_denoize -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
        | vsearch \
        --derep_fulllength - \
        --output $output_dir/zOTUs.fasta \
        --uc $output_dir/zOTUs.uc \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout
      fi
      
    fi # end of individual denoising

    ### OTU tables
    ### Cat dereplicated individual samples for making an OTU table
    cat dereplicated_sequences/*.fasta > $output_dir/Dereplicated_samples.fasta

    ## Prepare table with sequence abundance per sample
    seqkit seq --name $output_dir/Dereplicated_samples.fasta \
      | awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
      | sed 's/size=//; s/sample=//' \
      > tempdir/ASV_table_long.txt

    ## zOTU table creation
    printf "Making zOTU table ... \n"
    Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
      --derepuc      tempdir/Glob_derep.uc \
      --uc           "$output_dir"/zOTUs.uc \
      --asv          tempdir/ASV_table_long.txt \
      --rmsingletons FALSE \
      --fasta        $output_dir/zOTUs.fasta \
      --output       "$output_dir"/zOTU_table.temp 2>&1)
    echo $Rlog > tempdir/zOTU_table_creation.log 
    wait

    #format R-log file
    sed -i "s/;; /\n/g" $output_dir/zOTU_table_creation.log

    ### remove ";sample=.*;" and ";size=" from zOTUs.fasta file.
    if [[ -f $output_dir/zOTUs.fasta ]]; then
      sed -i 's/;sample=.*//' $output_dir/zOTUs.fasta
    fi

    # match zOTUs in a zOTU table with the IDs in the zOTUs.fasta file (may be different because of chimera removal step)
    checkerror=$(seqkit seq \
                -n $output_dir/zOTUs.fasta \
                > $output_dir/zOTUs_IDs.txt 2>&1)
    check_app_error
    checkerror=$(grep -f $output_dir/zOTUs_IDs.txt $output_dir/zOTU_table.temp \
                        > $output_dir/zOTU_table.txt 2>&1)
    check_app_error
    # remove intermediate files
    rm $output_dir/zOTUs_IDs.txt
    # add 1st row of the $input_table to the $output_dir/${input_table_base_name%%.txt}_lenFilt.temp
    header=$(head -n 1 $output_dir/zOTU_table.temp)
    sed -i "1i\\$header" "$output_dir/zOTU_table.txt"

    # remove intermediate files
    if [[ -f $output_dir/zOTU_table.temp ]]; then
        rm $output_dir/zOTU_table.temp
    fi

    ## Perform OTU clustering (if required, id < 1)
    if [[ $id_float != 1 ]]; then
      printf "\n Clustering zOTUs ... \n"

      ### Clustering
      checkerror=$(vsearch --cluster_size \
      $output_dir/zOTUs.fasta \
      $id \
      $simtype \
      $strands \
      $mask \
      $maxaccepts \
      $maxrejects \
      $cores \
      --centroids $output_dir/OTUs.fasta \
      --uc $output_dir/OTUs.uc \
      --fasta_width 0 \
      --sizein --sizeout 2>&1)
      check_app_error

      ## OTU table creation
      printf "Making OTU table ... \n"
      Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
        --derepuc      tempdir/Glob_derep.uc \
        --uc           "$output_dir"/OTUs.uc \
        --asv          tempdir/ASV_table_long.txt \
        --rmsingletons FALSE \
        --fasta        $output_dir/OTUs.fasta \
        --output       "$output_dir"/OTU_table.txt 2>&1)
      echo $Rlog > tempdir/OTU_table_creation.log 
      wait

      # Store output files in arrays for multiple runs
      output_feature_tables2+=("$output_dir/OTU_table.txt")
      output_fastas2+=("$output_dir/OTUs.fasta")
    fi # end of OTU clustering

    ### remove ";sample=.*;" and ";size=" from OTU.fasta files.
        # removing ";size=" because OTU table does not have "size" annotations; so the files would fit to LULU
    if [[ -f $output_dir/OTUs.fasta ]]; then
      sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
      sed -i 's/;size=.*//' $output_dir/OTUs.fasta
    fi

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    printf "\nCleaning up and compiling final stats files ... \n"

    mv $output_dir/Glob_derep.fasta tempdir/Glob_derep.fasta
    mv $output_dir/Dereplicated_samples.fasta tempdir/Dereplicated_samples.fasta

    #If input = FASTQ, then mkdir for converted fasta files
    if [[ $was_fastq == "true" ]]; then
        mkdir -p $output_dir/clustering_input_to_FASTA
        mv *.fasta $output_dir/clustering_input_to_FASTA
    fi

    #Delete decompressed files if original set of files were compressed
    if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
      delete_uncompressed=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $(NF-1)}';))
      rm *.$delete_uncompressed
    fi

    #Delete tempdirs
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir ]]; then
            rm -rf tempdir
        fi
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        if [[ -d tempdir_denoize ]]; then
            rm -rf tempdir_denoize
        fi
        if [[ -d tempdir_chimera ]]; then
            rm -rf tempdir_chimera
        fi
        if [[ -f $output_dir/zOTU_table_creation.log ]]; then
            rm -f $output_dir/zOTU_table_creation.log
        fi
        if [[ -f $output_dir/OTU_table_creation.log ]]; then
            rm -f $output_dir/OTU_table_creation.log
        fi
      else 
        #compress files in /tempdir
        if [[ -d tempdir ]]; then
            pigz tempdir/*
        fi
    fi

    #Make README.txt file
    # count features and sequences; outputs variables feature_count, nSeqs, nSample
    count_features "$output_dir/zOTU_table.txt"

    end=$(date +%s)
    runtime=$((end-start))

    printf "# Denoising was performed using UNOISE3.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'clustering_out' directory:
-----------------------------------
# zOTUs.fasta    = FASTA formated denoized sequences (zOTUs.fasta). Headers are renamed according to sha1 algorithm in vsearch.
# zOTU_table.txt = zOTU distribution table per sample (per input file in the working directory).
# zOTUs.uc       = uclust-like formatted clustering results for zOTUs

Number of zOTUs (zero-radius OTUs)    = $feature_count
Number of sequences in the zOTU table = $nSeqs
Number of samples in the zOTU table   = $nSample

Core command -> 
UNOISE: vsearch --cluster_unoise dereplicated_sequences.fasta $strands $minsize $unoise_alpha $simtype $mask $maxaccepts $maxrejects $cores --centroids zOTUs.fasta --uc zOTUs.uc \n\n" > $output_dir/README.txt

    ## If additional clustering was performed
    if [[ $id_float != 1 ]]; then
        size_otu=$(grep -c "^>" $output_dir/OTUs.fasta)
        count_features "$output_dir/OTU_table.txt"
        printf "Additional clustering of zOTUs at $id similarity threshold formed $size_otu OTUs.
    # OTUs.fasta    = FASTA formated representative OTU sequences. Headers are renamed according to sha1 algorithm in vsearch.
    # OTU_table.txt = OTU distribution table per sample (per input file in the working directory).
    # OTUs.uc       = uclust-like formatted clustering results for OTUs.
    
    Number of Features                       = $feature_count
    Number of sequences in the Feature table = $nSeqs
    Number of samples in the Feature table   = $nSample

    Core command -> 
    clustering: vsearch --cluster_size zOTUs.fasta $id $simtype $strands $mask $maxaccepts $maxrejects $cores --centroids OTUs.fasta --uc OTUs.uc \n\n" >> $output_dir/README.txt
    fi

    ## Chimera stats
    if [[ $chimerarm == "true" ]]; then
        printf "Chimera removal step eliminated $chimeras sequences\n" >> $output_dir/README.txt
    fi

    ## if input was fastq
    if [[ $was_fastq == "true" ]]; then
      printf "\nInput was fastq; converted those to fasta before clustering. 
      Converted fasta files in directory 'clustering_input_to_FASTA' \n" >> $output_dir/README.txt
    fi

    printf "\nIf samples are denoised individually rather by pooling all samples together, 
    reducing minsize to 4 is more reasonable for higher sensitivity.
    \n" >> $output_dir/README.txt

    printf "\n
    ##############################################
    ###Third-party applications for this process:
    #vsearch v2.23.0 for clustering
        #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
        #https://github.com/torognes/vsearch
    #GNU Parallel 20210422 for job parallelisation 
        #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
    ################################################" >> $output_dir/README.txt
    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done


# Validate and combine output files
if [[ $multiDir == "TRUE" ]]; then
    # Create an array of dereplicated_sequences directories (for merge_runs_unoise_wf.sh)
    derep_dirs=()
    for dir in $DIRS; do
        derep_dir="${dir}/dereplicated_sequences"
        derep_dirs+=("$derep_dir")
    done
    printf "%s\n" "${derep_dirs[@]}" > /input/multiRunDir/.derep_seqs_dirs


    # Check each file in the arrays
    valid_feature_tables=()
    valid_fastas=()
    valid_feature_tables2=()
    valid_fastas2=()
    
    for table in "${output_feature_tables[@]}"; do
        if [[ -f "$table" && -s "$table" ]]; then
            valid_feature_tables+=("$table")
        else
            printf "Warning: zOTU table not found or empty: $table\n" >&2
        fi
    done
    
    for fasta in "${output_fastas[@]}"; do
        if [[ -f "$fasta" && -s "$fasta" ]]; then
            valid_fastas+=("$fasta")
        else
            printf "Warning: zOTU FASTA file not found or empty: $fasta\n" >&2
        fi
    done

    # Check OTU clustering outputs if performed
    if [[ $id_float != 1 ]]; then
        for table in "${output_feature_tables2[@]}"; do
            if [[ -f "$table" && -s "$table" ]]; then
                valid_feature_tables2+=("$table")
            else
                printf "Warning: OTU table not found or empty: $table\n" >&2
            fi
        done
        
        for fasta in "${output_fastas2[@]}"; do
            if [[ -f "$fasta" && -s "$fasta" ]]; then
                valid_fastas2+=("$fasta")
            else
                printf "Warning: OTU FASTA file not found or empty: $fasta\n" >&2
            fi
        done
    fi
      
    output_feature_table=$(IFS=,; echo "${valid_feature_tables[*]}")
    output_fasta=$(IFS=,; echo "${valid_fastas[*]}")
    if [[ $id_float != 1 ]]; then
        output_feature_table2=$(IFS=,; echo "${valid_feature_tables2[*]}")
        output_fasta2=$(IFS=,; echo "${valid_fastas2[*]}")
    fi
else
    # Check single run output files
    if [[ ! -f "$output_feature_table" || ! -s "$output_feature_table" ]]; then
        printf "Error: zOTU table not found or empty: $output_feature_table\n" >&2
        exit 1
    fi
    if [[ ! -f "$output_fasta" || ! -s "$output_fasta" ]]; then
        printf "Error: zOTUs FASTA file not found or empty: $output_fasta\n" >&2
        exit 1
    fi
    # Check OTU clustering outputs if performed
    if [[ $id_float != 1 ]]; then
        if [[ ! -f "$output_feature_table2" || ! -s "$output_feature_table2" ]]; then
            printf "Error: OTU table not found or empty: $output_feature_table2\n" >&2
            exit 1
        fi
        if [[ ! -f "$output_fasta2" || ! -s "$output_fasta2" ]]; then
            printf "Error: OTU FASTA file not found or empty: $output_fasta2\n" >&2
            exit 1
        fi
    fi
fi

if [[ $multiDir == "TRUE" ]]; then
# write clustering parameters into file for the merge_runs_vsearch_wf.sh
    cat > /input/multiRunDir/.clustering_params << EOF
id="${id}"
id_float="${id_float}"
strands="${strands}"
chimerarm="${chimerarm}"
denoise_level="${denoise_level}"
unoise_alpha="${unoise_alpha}"
minsize="${minsize}"
cores="${cores}"
abskew="${abskew}"
simtype="${simtype}"
maxaccepts="${maxaccepts}"
maxrejects="${maxrejects}"
mask="${mask}"
EOF
fi

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"
if [[ $id_float != 1 ]]; then
    echo "output_feature_table2=$output_feature_table2"
    echo "output_fasta2=$output_fasta2"
fi

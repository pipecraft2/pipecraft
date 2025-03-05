#!/bin/bash

# Merge sequencing runs processed with UNOISE workflow if working with multuple runs in multiRunDir. 
 # UNOISE clustering for all fasta files in either qualFiltered_out or ITSx_out directories.
  # + apply 'curate otu table' when this is applied.
 # Samples with the same name across runs are not automatically merged together.

 # 1. clustering all samples from all runs.
 # 2. Split tables per run
 # 3. curate tables (tj + lenFilt) if 'curate otu table' is enabled
 # 4. merge zOTU tables (and also OTU tables if present)
################################################
###Third-party applications:
# vsearch v2.29.4
# seqkit v2.9.0
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')    
printf "# vsearch version: $vsearch_version\n"
printf "# seqkit version: $seqkit_version\n"

start_time=$(date)
start=$(date +%s)

# Source the clustering parameters from previous step
if [[ -f "/input/multiRunDir/.clustering_params" ]]; then
    source /input/multiRunDir/.clustering_params
    printf "\n UNOISE3 parameters:\n"
    printf "id: $id\n"
    printf "id_float: $id_float\n"
    printf "strands: $strands\n"
    printf "chimerarm: $chimerarm\n"
    printf "denoise_level: $denoise_level\n"
    printf "unoise_alpha: $unoise_alpha\n"
    printf "minsize: $minsize\n"
    printf "cores: $cores\n"
    printf "abskew: $abskew\n"
    printf "simtype: $simtype\n"
    printf "maxaccepts: $maxaccepts\n"
    printf "maxrejects: $maxrejects\n"
    printf "mask: $mask\n"

else
    echo "Error: Could not find UNOISE3 parameters file" >&2
    exit 1
fi

# Source the curate table parameters from previous step
if [[ -f "/input/multiRunDir/.curate_table_params" ]]; then
    source /input/multiRunDir/.curate_table_params
    printf "\nCurate table parameters:\n"
    printf "f_value: $f_value\n" 
    printf "p_value: $p_value\n"
    printf "collapseNoMismatch: $collapseNoMismatch\n"
    printf "max_length: $max_length\n"
    printf "min_length: $min_length\n"
    printf "max_length_seqkit: $max_length_seqkit\n"
    printf "min_length_seqkit: $min_length_seqkit\n"
    printf "min_length_num: $min_length_num\n"
    printf "max_length_num: $max_length_num\n"
    curate_otu_table="true"

    # Set filter_tag_jumps based on p_value and f_value
    if [[ -n "$p_value" && "$p_value" != "0" ]] && [[ -n "$f_value" && "$f_value" != "0" ]]; then
        filter_tag_jumps="true"
        printf "filter_tag_jumps: $filter_tag_jumps\n"
    else
        filter_tag_jumps="false"
        printf "filter_tag_jumps: $filter_tag_jumps\n"
    fi
else
    echo "CURATE OTU TABLE = FALSE"
    curate_otu_table="false"
fi

# source for functions
source /scripts/submodules/framework.functions.sh

# Excecute only if multiDir = true
if [[ ! -d "/input/multiRunDir" ]] && [[ $merge_runs != "true" ]]; then
    printf '%s\n' "ERROR]: multiRunDir not detected. Cannot merge sequencing runs since you are working with a single sequencing run. 
    >DONE." >&2
    end_process
elif [[ $merge_runs == "true" ]]; then
    printf "\nStarting to merge sequencing runs...\n"
    cd /input/multiRunDir
    #output dir
    output_dir=$"/input/multiRunDir/merged_runs"
    export output_dir
    # remove output dir if it already exists
    if [[ -d "$output_dir" ]]; then
        rm -rf $output_dir
    fi
    # create new output dir
    mkdir -p $output_dir
    echo "output dir: $output_dir"
else
    printf '%s\n' "ERROR]: Merge sequencing runs is not enabled. Exiting.\n" >&2
    end_process
fi

###################################################
# 1. UNOISE clustering all samples from all runs. #
###################################################
 # using fasta files qualFiltered_out (or ITSx_out if ITSx = true); 
 # dereplicated_sequences dirs are exported to file /input/multiRunDir/.derep_seqs_dirs by clustering_unoise.sh (this dir capture is solved in clustering_unoise.sh)

# Check if .derep_seqs_dirs exists
if [[ ! -f "/input/multiRunDir/.derep_seqs_dirs" ]]; then
    printf '%s\n' "ERROR: /input/multiRunDir/.derep_seqs_dirs not found. This file should contain paths to dereplicated sequence directories.
    >DONE." >&2
    end_process
fi

# create dereplicated_sequences directory
if [[ -d "$output_dir/dereplicated_sequences" ]]; then
    rm -rf $output_dir/dereplicated_sequences
fi
mkdir -p $output_dir/dereplicated_sequences

# Create a directory for run sample lists
if [[ -d "$output_dir/run_sample_lists" ]]; then
    rm -rf $output_dir/run_sample_lists
fi
mkdir -p $output_dir/run_sample_lists

# Process each directory and track which samples belong to which run
while IFS= read -r dir_path; do
    if [[ -d "$dir_path" ]]; then
        printf "Copying fasta files from: %s\n" "$dir_path"
        
        # Extract run name from the directory path (gets the first directory component)
        run=$(echo "$dir_path" | cut -d'/' -f1)
        echo "Run name: $run"
        
        # Copy the fasta files
        cp "$dir_path"/*.fasta "$output_dir/dereplicated_sequences"
        
        # Create a list of sample names for this run
        ls "$dir_path"/*.fasta | sed 's/.*\///' | sed 's/\.fasta$//' > "$output_dir/run_sample_lists/${run}_samples.txt"
    else
        printf "Warning: Directory not found: %s\n" "$dir_path"
    fi
done < "/input/multiRunDir/.derep_seqs_dirs"

### Global dereplication
if [[ -d $output_dir/tempdir ]]; then
    rm -rf $output_dir/tempdir
fi
mkdir -p $output_dir/tempdir

printf "Dereplicating globally ... \n"
find $output_dir/dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
| vsearch \
--derep_fulllength - \
--uc $output_dir/tempdir/Glob_derep.uc \
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
    find $output_dir/dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "denoise_and_chim {}"

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
cat $output_dir/dereplicated_sequences/*.fasta > $output_dir/Dereplicated_samples.fasta

## Prepare table with sequence abundance per sample
seqkit seq --name $output_dir/Dereplicated_samples.fasta \
    | awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
    | sed 's/size=//; s/sample=//' \
    > $output_dir/tempdir/ASV_table_long.txt

## zOTU table creation
printf "Making zOTU table ... \n"
Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
    --derepuc      $output_dir/tempdir/Glob_derep.uc \
    --uc           $output_dir/zOTUs.uc \
    --asv          $output_dir/tempdir/ASV_table_long.txt \
    --rmsingletons FALSE \
    --fasta        $output_dir/zOTUs.fasta \
    --output       $output_dir/zOTU_table.txt 2>&1)
echo $Rlog > $output_dir/tempdir/zOTU_table_creation.log 
wait

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
    --derepuc      $output_dir/tempdir/Glob_derep.uc \
    --uc           $output_dir/OTUs.uc \
    --asv          $output_dir/tempdir/ASV_table_long.txt \
    --rmsingletons FALSE \
    --fasta        $output_dir/OTUs.fasta \
    --output       $output_dir/OTU_table.txt 2>&1)
    echo $Rlog > $output_dir/tempdir/OTU_table_creation.log 
    wait
fi # end of OTU clustering

### remove ";sample=.*;" and ";size=" from OTU.fasta files.
    # removing ";size=" because OTU table does not have "size" annotations; so the files would fit to LULU
if [[ -f $output_dir/OTUs.fasta ]]; then
    sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta
fi
if [[ -f $output_dir/zOTUs.fasta ]]; then
    sed -i 's/;sample=.*;/;/' $output_dir/zOTUs.fasta
    sed -i 's/;size=.*//' $output_dir/zOTUs.fasta
fi

###########################
# 2. Split tables per run #
###########################
# outputs per run OTU/zOTU table and fasta file in split_tables directory
# Global OTU/zOTU table must have 'Sequence' as 2nd column
printf "Splitting zOTU/OTU tables per run...\n"

# Get run names from multiRunDir
cd /input/multiRunDir
RUNS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")

# Create directory for split tables
if [[ -d "$output_dir/split_tables" ]]; then
    rm -rf $output_dir/split_tables
fi
mkdir -p $output_dir/split_tables

# Unset arrays if they exist
unset output_feature_tables
unset output_fastas
unset output_zfeature_tables
unset output_zfastas

# Create output arrays
declare -a output_zfeature_tables
declare -a output_zfastas
declare -a output_feature_tables
declare -a output_fastas

# Function to split table by run
split_table_by_run() {
    local input_table=$1
    local input_fasta=$2
    local table_prefix=$3
    local fasta_prefix=$4
    local run=$5

    # Get sample names for this run from the run_sample_lists directory
    sample_list=$(printf "OTU\nSequence\n"; cat "$output_dir/run_sample_lists/${run}_samples.txt")
    # Get column numbers and join them with commas (including OTU column)
    cols=$(head -n1 $input_table | tr '\t' '\n' | nl -v0 | grep -f <(echo "$sample_list") | cut -f1 | awk '{print $1 + 1}' | paste -sd',')
    # Extract the columns using cut with the column numbers
    cut -f"$cols" $input_table > $output_dir/split_tables/${table_prefix}_${run}.temp
    
    # Filter zero-abundance rows (from col 3-...) and export both table and FASTA
    table_file=$output_dir/split_tables/${table_prefix}_${run}.temp
    awk -v outdir="$output_dir/split_tables" -v run="$run" -v table_prefix="$table_prefix" -v fasta_prefix="$fasta_prefix" '
        NR==1{print $0 > outdir"/"table_prefix"_"run".txt"}
        NR>1{
            sum=0
            for(i=3;i<=NF;i++) sum+=$i
            if(sum>0) {
                print $0 > outdir"/"table_prefix"_"run".txt"
                print ">"$1"\n"$2 > outdir"/"fasta_prefix"_"run".fasta"
            }
    }' "$table_file"
    rm $table_file
}

# For each run, extract samples belonging to that run from the OTU/zOTU tables
for run in $RUNS; do
    printf "Processing tables for $run...\n"
    
    # Process zOTU table if it exists
    if [[ -f "$output_dir/zOTU_table.txt" ]]; then
        printf "Getting zOTU table for $run...\n"
        split_table_by_run "$output_dir/zOTU_table.txt" "$output_dir/zOTUs.fasta" "zOTU_table" "zOTUs" "$run"
        
        # Add current run's outputs to arrays
        output_zfeature_tables+=("$output_dir/split_tables/zOTU_table_${run}.txt")
        output_zfastas+=("$output_dir/split_tables/zOTUs_${run}.fasta")
    fi
    # Process OTU table if it exists
    if [[ -f "$output_dir/OTU_table.txt" ]]; then
        printf "Getting OTU table for $run...\n"
        split_table_by_run "$output_dir/OTU_table.txt" "$output_dir/OTUs.fasta" "OTU_table" "OTUs" "$run"
        
        # Add current run's outputs to arrays
        output_feature_tables+=("$output_dir/split_tables/OTU_table_${run}.txt")
        output_fastas+=("$output_dir/split_tables/OTUs_${run}.fasta")
    fi
done

# Print the output arrays (inputs for the next steps)
if [[ ${#output_zfeature_tables[@]} -gt 0 ]]; then
    echo "Split zOTU tables per run: output_zfeature_tables = ${output_zfeature_tables[*]}"
    echo "zOTU Fasta files per run: output_zfastas = ${output_zfastas[*]}"
fi
if [[ ${#output_feature_tables[@]} -gt 0 ]]; then
    echo "Split OTU tables per run: output_feature_tables = ${output_feature_tables[*]}"
    echo "OTU Fasta files per run: output_fastas = ${output_fastas[*]}"
fi

####################################################################
# 3. curate tables (tj + lenFilt) if 'curate otu table' is enabled #
####################################################################
# make curated output dir if filter_tag_jumps is true and/or min_length_num or max_length_num is not 0 or empty
if [[ $curate_otu_table == "true" && \
    ( $filter_tag_jumps == "true" || \
      ( -n "$min_length_num" && "$min_length_num" != "0" ) || \
      ( -n "$max_length_num" && "$max_length_num" != "0" ) ) ]]; then
    # create subdirectory for curated resultswhile preserving original output_dir
    curated_dir="${output_dir}/split_tables/curated"
    if [[ -d "$curated_dir" ]]; then
        rm -rf "$curated_dir"
    fi
    mkdir -p "$curated_dir"
    export curated_dir

else # no curation is performed
    curated_dir="$output_dir/split_tables"
fi

################
### UNCROSS2 ###
################
# Process both OTU and zOTU tables if they exist
for table_type in "zOTU" "OTU"; do
    # Skip this type if no corresponding tables exist
    if [[ $table_type == "OTU" && ${#output_feature_tables[@]} -eq 0 ]] || \
       [[ $table_type == "zOTU" && ${#output_zfeature_tables[@]} -eq 0 ]]; then
        continue
    fi

    # Set the appropriate array based on table type
    if [[ $table_type == "OTU" ]]; then
        table_files=("${output_feature_tables[@]}")
        fasta_files=("${output_fastas[@]}")
    else
        table_files=("${output_zfeature_tables[@]}")
        fasta_files=("${output_zfastas[@]}")
    fi

    for ((i=0; i<${#table_files[@]}; i++)); do
        table_file="${table_files[i]}"
        fasta_file="${fasta_files[i]}"
        
        # get run name from table file name
        run=$(basename "$table_file" | sed "s/${table_type}_table_\(.*\)\.txt/\1/")
        echo "Processing $table_type table for run = $run"

        ### Process samples with UNCROSS2 (tag-jumps filtering) in R
        if [[ $filter_tag_jumps == "true" ]]; then
            table_file_basename=$(basename $table_file)    

            printf "# Running tag-jumps filtering (UNCROSS2) for $table_file_basename\n "
            # Filter primary feature table
            Rlog=$(Rscript /scripts/submodules/tag_jump_removal.R $table_file $f_value $p_value $fasta_file $curated_dir 2>&1)
            # Check if R script executed successfully
            if [ $? -ne 0 ]; then
                log_error "tag-jumps filtering R script failed with the following error:
                $Rlog
                Please check the parameters and input file.
                >Quitting"
                end_process
            fi
            echo "$Rlog" > "$curated_dir/tag-jumps_filt.log"

            # format R-log file
            sed -i "s/;; /\n/g" $curated_dir/tag-jumps_filt.log 
        
            # Check if output files were created
            if [ -z "$(find "$curated_dir" -name "*_TagJumpFilt.txt")" ]; then
                log_error "tag-jumps filtering process did not generate the expected output file.
                Please check the log file at $curated_dir/tag-jumps_filt.log
                >Quitting"
                end_process
            fi
            printf "tag-jumps filtering completed \n\n"

            # cp fasta file to curated_dir (merged_runs/split_tables/curated)
            cp "$fasta_file" "$curated_dir"

            # Update output variables if only tag-jumps filtering is performed
            if [[ -z $min_length_num || $min_length_num == "0" ]] && [[ -z $max_length_num || $max_length_num == "0" ]]; then
                printf "Only tag-jumps filtering is performed\n"
                if [[ $table_type == "OTU" ]]; then
                    # Update OTU arrays
                    output_feature_tables[$i]=$curated_dir/${table_file_basename%%.txt}_TagJumpFilt.txt
                    output_fastas[$i]=$curated_dir/${table_type}s_${run}.fasta
                else
                    # Update zOTU arrays
                    output_zfeature_tables[$i]=$curated_dir/${table_file_basename%%.txt}_TagJumpFilt.txt
                    output_zfastas[$i]=$curated_dir/${table_type}s_${run}.fasta
                fi

                # Print updated arrays for debugging
                if [[ $table_type == "OTU" ]]; then
                    echo "Updated OTU arrays:"
                    echo "output_feature_tables[$i]=${output_feature_tables[$i]}"
                    echo "output_fastas[$i]=${output_fastas[$i]}"
                else
                    echo "Updated zOTU arrays:"
                    echo "output_zfeature_tables[$i]=${output_zfeature_tables[$i]}"
                    echo "output_zfastas[$i]=${output_zfastas[$i]}"
                fi

            else 
                printf "Length filtering is also performed\n"
                # Set inputs for length filtering
                input_table=$curated_dir/${table_file_basename%%.txt}_TagJumpFilt.txt
                input_fasta=$curated_dir/${table_type}s_${run}.fasta
            fi

        ### skip tag-jumps filtering ###
        elif [[ $filter_tag_jumps == "false" ]]; then
            printf "# Skipping tag-jumps filtering\n"
            input_table=$table_file
            input_fasta=$fasta_file
        fi
   
        ########################
        ### length filtering ###
        ########################
        if  [[ $min_length_num != "0" && -n $min_length_num ]] || \
            [[ $max_length_num != "0" && -n $max_length_num ]]; then
            printf "Filtering by length, min_length = $min_length_num, max_length = $max_length_num. \n"
            # get basenames for correct naming of output files
            input_table_basename=$(basename $input_table)
            fasta_basename=$(basename $input_fasta)
            # count input ASVs
            ASVs_count=$(grep -c "^>" $input_fasta)
            echo "Feature count = $ASVs_count"
            # filter by length
            checkerror=$(seqkit seq -w 0 -g \
                        $min_length_seqkit \
                        $max_length_seqkit \
                        $input_fasta \
                        > $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta 2>&1)
            check_app_error

            # count length filtered ASVs and proceed with the rest of the steps
            ASVs_lenFilt=$(grep -c "^>" $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta)
            echo "length filtered Feature (ASV/OTU) count = $ASVs_lenFilt"
            if [[ $ASVs_lenFilt == 0 ]]; then
                ASVs_lenFilt_result=$"All Features (ASVs/OTUs) were filtered out based on the length filter
                (min_length $min_length_num bp and max_length $max_length_num bp).
                No new files generated.
                Input table was $input_table_basename and input fasta was $fasta_basename with $ASVs_count sequences"
                echo -e "$ASVs_lenFilt_result"
                rm $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta
                # Set output variables to input files since no filtering occurred
                if [[ $table_type == "OTU" ]]; then
                    output_feature_table=$input_table
                    output_fasta=$input_fasta
                else
                    output_zfeature_table=$input_table
                    output_zfasta=$input_fasta
                fi
                
            elif [[ $ASVs_lenFilt == $ASVs_count ]]; then
                ASVs_lenFilt_result=$"None of the Features (ASVs/OTUs) were filtered out based on the length filter
                (min_length $min_length_num bp and max_length $max_length_num bp).
                No new files generated.
                Input table was $input_table_basename and input fasta was $fasta_basename with $ASVs_count sequences"
                echo -e "$ASVs_lenFilt_result"
                export ASVs_lenFilt_result
                rm $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta
                # Set output variables to input files since no filtering occurred
                if [[ $table_type == "OTU" ]]; then
                    output_feature_table=$input_table
                    output_fasta=$input_fasta
                else
                    output_zfeature_table=$input_table
                    output_zfasta=$input_fasta
                fi
            else          
                # filter the table
                checkerror=$(seqkit seq \
                            -n $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta \
                            > $curated_dir/${fasta_basename%%.fasta}_IDs.txt 2>&1)
                check_app_error
                checkerror=$(grep -f $curated_dir/${fasta_basename%%.fasta}_IDs.txt $input_table \
                                    > $curated_dir/${input_table_basename%%.txt}_lenFilt.temp 2>&1)
                check_app_error
                # remove intermediate files
                rm $curated_dir/${fasta_basename%%.fasta}_IDs.txt
                # add 1st row of the $input_table to the $output_dir/${input_table_basename%%.txt}_lenFilt.temp
                header=$(head -n 1 $input_table)
                sed -i "1i\\$header" "$curated_dir/${input_table_basename%%.txt}_lenFilt.temp"

                # Remove samples (columns) with 0 abundance; does not remove 0 rows, but there cannot be 0 rows
                checkerror=$(awk '
                BEGIN {
                    FS = OFS = "\t"
                }
                NR == 1 {
                    # Store the header and identify the "Sequence" column
                    for (i = 1; i <= NF; i++) {
                        headers[i] = $i
                        if ($i == "Sequence") {
                            sequence_col = i
                        }
                    }
                    next
                }
                {
                    # Sum each column, excluding the first column and the "Sequence" column
                    for (i = 2; i <= NF; i++) {
                        if (i != sequence_col) {
                            sum[i] += $i
                        }
                    }
                    # Store the row data
                    for (i = 1; i <= NF; i++) {
                        data[NR, i] = $i
                    }
                }
                END {
                    # Print the header with non-zero columns
                    for (i = 1; i <= NF; i++) {
                        if (i == 1 || i == sequence_col || sum[i] != 0) {
                            printf "%s%s", headers[i], (i == NF ? "\n" : OFS)
                        }
                    }
                    # Print the rows excluding columns with zero sum
                    for (j = 2; j <= NR; j++) {
                        for (i = 1; i <= NF; i++) {
                            if (i == 1 || i == sequence_col || sum[i] != 0) {
                                printf "%s%s", data[j, i], (i == NF ? "\n" : OFS)
                            }
                        }
                    }
                }
                ' "$curated_dir/${input_table_basename%%.txt}_lenFilt.temp" \
                > "$curated_dir/${input_table_basename%%.txt}_lenFilt.txt" 2>&1)
                check_app_error
                # remove intermediate files
                if [[ -f $curated_dir/${input_table_basename%%.txt}_lenFilt.temp ]]; then
                    rm $curated_dir/${input_table_basename%%.txt}_lenFilt.temp
                fi
                # for the report
                count_features "$curated_dir/${input_table_basename%%.txt}_lenFilt.txt"
                ASVs_lenFilt_result=$"Features (ASVs/OTUs) after length filtering = $ASVs_lenFilt.

- ${input_table_basename%%.txt}_lenFilt.txt = Feature table after length filtering.
- ${fasta_basename%%.fasta}_lenFilt.fasta = Representative sequences file after length filtering

Number of Features                       = $feature_count
Number of sequences in the Feature table = $nSeqs
Number of samples in the Feature table   = $nSample"
                echo -e "$ASVs_lenFilt_result"
                # Set output variables to filtered files
                if [[ $table_type == "OTU" ]]; then
                    output_feature_table="$curated_dir/${input_table_basename%%.txt}_lenFilt.txt"
                    output_fasta="$curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta"
                else
                    output_zfeature_table="$curated_dir/${input_table_basename%%.txt}_lenFilt.txt"
                    output_zfasta="$curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta"
                fi
            fi
        fi

        # Update the appropriate arrays with the final output files for this run
        if [[ $table_type == "OTU" ]]; then
            output_feature_tables[$i]=$output_feature_table
            output_fastas[$i]=$output_fasta
        else
            output_zfeature_tables[$i]=$output_zfeature_table
            output_zfastas[$i]=$output_zfasta
        fi

    done
    
    # Print final arrays for this table type
    if [[ $table_type == "OTU" ]]; then
        echo "OTU tables to be merged = ${output_feature_tables[*]}"
        echo "Corresponding OTU fasta files = ${output_fastas[*]}"
    else
        echo "zOTU tables to be merged = ${output_zfeature_tables[*]}"
        echo "Corresponding zOTU fasta files = ${output_zfastas[*]}"
    fi
done

################################
# 4. merge zOTU and OTU tables #
################################
# Merge zOTU tables if they exist
if [[ ${#output_zfeature_tables[@]} -gt 0 ]]; then
    printf "Merging zOTU tables...\n"
    output_feature_table=$(IFS=,; echo "${output_zfeature_tables[*]}")
    output_fasta=$(IFS=,; echo "${output_zfastas[*]}")
    export output_feature_table
    export output_fasta

    # Merge tables in R 
    Rlog=$(Rscript /scripts/submodules/mergeOtuTables.R $output_feature_table $output_dir zOTU_table.txt 2>&1)
    # Check if R script executed successfully
    if [ $? -ne 0 ]; then
        log_error "mergeOtuTables R script failed with the following error:
        $Rlog
        Please check the parameters and input file.
        >Quitting"
        end_process
    fi
    echo "$Rlog" > "$output_dir/mergedOtuTables.log"

    # Merge FASTA files
    printf "Merging zOTU FASTA files...\n"
    cat "${output_zfastas[@]}" > "$output_dir/merged_zOTUs.fasta"
    # Remove duplicate sequences from merged FASTA file
    printf "Removing duplicate sequences from merged zOTU FASTA file...\n"
    checkerror=$(seqkit rmdup -w 0 -s "$output_dir/merged_zOTUs.fasta" -o "$output_dir/merged_zOTUs_dedup.fasta" 2>&1)
    check_app_error
    mv "$output_dir/merged_zOTUs_dedup.fasta" "$output_dir/zOTUs.fasta" && rm "$output_dir/merged_zOTUs.fasta"
fi

# Merge OTU tables if they exist
if [[ ${#output_feature_tables[@]} -gt 0 ]]; then
    printf "Merging OTU tables...\n"
    output_feature_table=$(IFS=,; echo "${output_feature_tables[*]}")
    output_fasta=$(IFS=,; echo "${output_fastas[*]}")
    export output_feature_table
    export output_fasta

    # Merge tables in R 
    Rlog=$(Rscript /scripts/submodules/mergeOtuTables.R $output_feature_table $output_dir OTU_table.txt 2>&1)
    # Check if R script executed successfully
    if [ $? -ne 0 ]; then
        log_error "mergeOtuTables R script failed with the following error:
        $Rlog
        Please check the parameters and input file.
        >Quitting"
        end_process
    fi
    echo "$Rlog" > "$output_dir/mergeOtuTables.log"

    # Merge FASTA files
    printf "Merging OTU FASTA files...\n"
    cat "${output_fastas[@]}" > "$output_dir/merged_OTUs.fasta"
    # Remove duplicate sequences from merged FASTA file
    printf "Removing duplicate sequences from merged OTU FASTA file...\n"
    checkerror=$(seqkit rmdup -w 0 -s "$output_dir/merged_OTUs.fasta" -o "$output_dir/merged_OTUs_dedup.fasta" 2>&1)
    check_app_error
    mv "$output_dir/merged_OTUs_dedup.fasta" "$output_dir/OTUs.fasta" && rm "$output_dir/merged_OTUs.fasta"
fi

### CLEAN UP ###
# Remove all files and folders except OTUs.fasta, OTU_table.txt, zOTUs.fasta and zOTU_table.txt
printf "Cleaning up temporary files...\n"
if [[ $debugger != "true" ]]; then
    find "$output_dir" -mindepth 1 -not \( -name "OTUs.fasta" -o -name "OTU_table.txt" -o -name "zOTUs.fasta" -o -name "zOTU_table.txt" \) -exec rm -rf {} +

    if [[ -f "/input/multiRunDir/.curate_table_params" ]]; then
        rm /input/multiRunDir/.curate_table_params
    fi
    if [[ -f "/input/multiRunDir/.clustering_params" ]]; then
        rm /input/multiRunDir/.clustering_params
    fi
    if [[ -f "/input/multiRunDir/.derep_seqs_dirs" ]]; then
        rm /input/multiRunDir/.derep_seqs_dirs
    fi
fi 

##########################################
### Make README.txt file (merged_runs) ###
##########################################
# Count features for both tables if they exist
if [[ -f "$output_dir/zOTU_table.txt" ]]; then
    count_features "$output_dir/zOTU_table.txt"
    zfeature_count=$feature_count
    znSeqs=$nSeqs
    znSample=$nSample
fi

if [[ -f "$output_dir/OTU_table.txt" ]]; then
    count_features "$output_dir/OTU_table.txt"
    ofeature_count=$feature_count
    onSeqs=$nSeqs
    onSample=$nSample
fi

end=$(date +%s)
runtime=$((end-start))

printf "# Merged sequencing runs: 

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

# 1. clustering all samples from all sequencing runs: 
        similarity_threshold = $id
        strands = $strands
        denoise_level = $denoise_level
        unoise_alpha = $unoise_alpha
        minsize = $minsize
        similarity_type = $simtype
        maxaccepts = $maxaccepts
        mask = $mask
# 2. Split OTU/zOTU tables per run\n" > $output_dir/README.txt

if [[ $curate_otu_table == "true" ]]; then
    printf "# 3. curate OTU/zOTU tables per run:
    f_value = $f_value
    p_value = $p_value
    $min_length
    $max_length
# 4. merge (curated) OTU/zOTU tables\n" >> $output_dir/README.txt
else
    printf "# 3. merge OTU/zOTU tables\n" >> $output_dir/README.txt
fi

printf "\nOutput files:\n------------" >> $output_dir/README.txt

if [[ -f "$output_dir/zOTU_table.txt" ]]; then
    printf "\n# zOTU_table.txt = merged zOTU abundance table
# zOTUs.fasta    = merged zOTU sequences

Number of zOTUs                       = $zfeature_count
Number of sequences in zOTU table     = $znSeqs
Number of samples in zOTU table       = $znSample" >> $output_dir/README.txt
fi

if [[ -f "$output_dir/OTU_table.txt" ]]; then
    printf "\n\n# OTU_table.txt = merged OTU abundance table
# OTUs.fasta    = merged OTU sequences

Number of OTUs                       = $ofeature_count
Number of sequences in OTU table     = $onSeqs
Number of samples in OTU table       = $onSample" >> $output_dir/README.txt
fi

printf "\n
#############################################
###Third-party applications for this process:
# vsearch (version $vsearch_version)
    #citation: Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584
    #https://github.com/torognes/vsearch
# seqkit (version $seqkit_version)
    #citation: Shen W, Le S, Li Y, Hu F, Li Z, et al. (2016) SeqKit: A Multi-Format Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://github.com/shenwei356/seqkit
#############################################" >> $output_dir/README.txt

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

# variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"

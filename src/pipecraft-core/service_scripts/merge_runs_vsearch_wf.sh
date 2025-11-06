#!/bin/bash

# Merge sequencing runs processed with vsearch OTUs workflow if working with multuple runs in multiRunDir. 
 # vsearch clustering for all fasta files in either chimeraFiltered_out or ITSx_out directories.
  # + apply 'curate otu table' when this is applied.
 # Samples with the same name across runs are not automatically merged together.

 # 1. clustering all samples from all runs.
 # 2. Split tables per run
 # 3. curate tables (tj + lenFilt) if 'curate otu table' is enabled
 # 4. merge OTU tables
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
    printf "\nClustering parameters:\n"
    printf "similarity_threshold: $id\n" 
    printf "OTU_type: $otutype\n"
    printf "strands: $strands\n"
    printf "remove_singletons: $remove_singletons\n"
    printf "sequence_sorting: $seqsort\n"
    printf "similarity_type: $simtype\n"
    printf "centroid_in: $centroid_in\n"
    printf "maxaccepts: $maxaccepts\n"
    printf "mask: $mask\n"
    printf "cores: $cores\n"
else
    echo "Error: Could not find clustering parameters file" >&2
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

############################################
# 1. clustering all samples from all runs. #
############################################
 # using fasta files chimeraFiltered_out (or ITSx_out if ITSx = true); 
 # dereplicated_sequences dirs are exported to file /input/multiRunDir/.derep_seqs_dirs by clustering_vsearch.sh

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
    --output $output_dir/Glob_derep.fasta \
    --uc $output_dir/tempdir/Glob_derep.uc \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout

### Clustering
printf "Clustering ... \n"
checkerror=$(vsearch $seqsort \
    $output_dir/Glob_derep.fasta \
    $id \
    $simtype \
    $strands \
    $mask \
    $centroid_in \
    $maxaccepts \
    $cores \
    $otutype $output_dir/OTUs.temp.fasta \
    --uc $output_dir/OTUs.uc \
    --fasta_width 0 \
    --sizein --sizeout 2>&1)
check_app_error


### Merge dereplicated FASTA files (for making OTU table) in chunks
batch_size=1000
batch_num=0
# Find all .fasta files and merge them in groups
find "$output_dir/dereplicated_sequences" -type f -name "*.fasta" | \
xargs -n $batch_size | while read batch; do
    batch_num=$((batch_num + 1))
    echo "Merging batch $batch_num..."
    echo "$batch" | xargs cat > "$output_dir/tempdir/Dereplicated_batch_${batch_num}.fasta"
done
# Merge all batch files into the final combined file
cat "$output_dir"/tempdir/Dereplicated_batch_*.fasta > "$output_dir"/tempdir/Dereplicated_samples.fasta
# clean up temporary batch files
rm "$output_dir"/tempdir/Dereplicated_batch_*.fasta

## Prepare table with sequence abundance per sample
checkerror=$(seqkit seq --name $output_dir/tempdir/Dereplicated_samples.fasta \
| awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
| sed 's/size=//; s/sample=//' \
> $output_dir/tempdir/ASV_table_long.txt 2>&1)
check_app_error

### OTU table creation
printf "Making OTU table ... \n"
Rlog=$(Rscript /scripts/submodules/make_OTU_long_table.R \
--derepuc      "$output_dir"/tempdir/Glob_derep.uc \
--uc           "$output_dir"/OTUs.uc \
--asv          "$output_dir"/tempdir/ASV_table_long.txt \
--rmsingletons $remove_singletons \
--output       "$output_dir"/OTU_table_long.txt 2>&1)
echo $Rlog > "$output_dir"/tempdir/OTU_table_creation.log 
wait

### Discard singleton OTUs
if [[ $remove_singletons == "true"  ]]; then
    printf "Discarding singletons from fasta file... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --sizeout --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error
    # remove ";sample=.*;" from OTU.fasta file.
    sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
    # removing ";size=" because OTU table does not have "size" annotations; so the files would fit to LULU
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta 
    mv $output_dir/OTUs.temp.fasta "$output_dir"/tempdir/
else
    sed -e 's/;sample=.*;/;/' $output_dir/OTUs.temp.fasta > $output_dir/OTUs.fasta
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta
    mv $output_dir/OTUs.temp.fasta "$output_dir"/tempdir/
fi

###########################
# 2. Split tables per run #
###########################
# outputs per run OTU table and fasta file in split_tables directory
printf "Splitting OTU table per run...\n"

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

# Create output arrays
declare -a output_feature_tables
declare -a output_fastas

# For each run, extract samples belonging to that run from the long format OTU table
for run in $RUNS; do
    printf "Getting OTU table for $run...\n"
    
    # Filter long format table for samples in this run and remove zero-abundance rows
    awk -v run="$run" -v outdir="$output_dir/split_tables" -v sample_file="$output_dir/run_sample_lists/${run}_samples.txt" '
        BEGIN {
            # Read sample list for this run
            while ((getline < sample_file) > 0) {
                samples[$0] = 1
            }
            close(sample_file)
        }
        NR==1 {
            # Print header
            print $0 > outdir"/OTU_table_long_"run".txt"
            next
        }
        {
            # Filter by SampleID (first column) and non-zero abundance (third column)
            if ($1 in samples && $3 > 0) {
                print $0 > outdir"/OTU_table_long_"run".txt"
                otus[$2] = 1  # Store unique OTU names
            }
        }
        END {
            # Write unique OTU list for FASTA extraction
            for (otu in otus) {
                print otu > outdir"/OTU_list_"run".txt"
            }
        }
    ' "$output_dir/OTU_table_long.txt"
    
    # Extract OTU sequences from main FASTA file using seqkit
    if [[ -f "$output_dir/split_tables/OTU_list_${run}.txt" ]]; then
        seqkit grep -f "$output_dir/split_tables/OTU_list_${run}.txt" -w 0 \
            "$output_dir/OTUs.fasta" -o "$output_dir/split_tables/OTUs_${run}.fasta"
        rm "$output_dir/split_tables/OTU_list_${run}.txt"
    else
        # No OTUs for this run, create empty FASTA
        touch "$output_dir/split_tables/OTUs_${run}.fasta"
    fi

    # Initialize output variables
    output_feature_table=$output_dir/split_tables/OTU_table_long_${run}.txt
    output_fasta=$output_dir/split_tables/OTUs_${run}.fasta

    # Add current run's outputs to arrays (for cases when no curation is performed)
    output_feature_tables+=("$output_feature_table")
    output_fastas+=("$output_fasta")
done

# Print the output arrays (inputs for the next steps)
echo "Split tables per run: output_feature_tables = ${output_feature_tables[*]}"
echo "Fasta files per run: output_fastas = ${output_fastas[*]}"

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
for table_file in /input/multiRunDir/merged_runs/split_tables/OTU_table_long_*.txt; do

    # get run name from table file name
    run=$(basename "$table_file" | sed 's/OTU_table_long_\(.*\)\.txt/\1/')
    echo "run = $run"
    
    # table_file corresponding fasta file
    fasta_file=$output_dir/split_tables/OTUs_${run}.fasta

    ### Process samples with UNCROSS2 (tag-jumps filtering) in R
    if [[ $filter_tag_jumps == "true" ]]; then
        table_file_basename=$(basename $table_file)    

        printf "# Running tag-jumps filtering (UNCROSS2) for $table_file_basename\n "
        # Filter primary feature table
        Rlog=$(Rscript /scripts/submodules/tag_jump_removal_long.R $table_file $f_value $p_value $curated_dir 2>&1)
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
        table_dir=$(dirname "$table_file")
        cp "$table_dir/OTUs_${run}.fasta" "$curated_dir"

        # Update output variables if only tag-jumps filtering is performed
        if [[ -z $min_length_num || $min_length_num == "0" ]] && [[ -z $max_length_num || $max_length_num == "0" ]]; then
            printf "Only tag-jumps filtering is performed\n"
            output_feature_table=$curated_dir/${table_file_basename%%.txt}_TagJumpFilt.txt
            output_fasta=$curated_dir/OTUs_${run}.fasta
            export output_feature_table
            export output_fasta

            printf "Table for merging = $output_feature_table\n"
            printf "Fasta file for merging = $output_fasta\n"
        else 
            printf "Length filtering is also performed\n"
            # Set inputs for length filtering
            input_table=$curated_dir/${table_file_basename%%.txt}_TagJumpFilt.txt
            input_fasta=$curated_dir/OTUs_${run}.fasta
            export input_table
            export input_fasta

            printf "Table for next step = $input_table\n"
            printf "Fasta file for next step = $input_fasta\n"
        fi

    ### skip tag-jumps filtering ###
    elif [[ $filter_tag_jumps == "false" ]]; then
        printf "# Skipping tag-jumps filtering\n"
        input_table=$table_file
        export input_table
        fasta_file=$output_dir/split_tables/OTUs_${run}.fasta
        export fasta_file
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
        echo "Feature (ASVs/OTUs) count = $ASVs_count"
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
            output_feature_table=$input_table
            output_fasta=$input_fasta
            
        elif [[ $ASVs_lenFilt == $ASVs_count ]]; then
            ASVs_lenFilt_result=$"None of the Features (ASVs/OTUs) were filtered out based on the length filter
            (min_length $min_length_num bp and max_length $max_length_num bp).
            No new files generated.
            Input table was $input_table_basename and input fasta was $fasta_basename with $ASVs_count sequences"
            echo -e "$ASVs_lenFilt_result"
            export ASVs_lenFilt_result
            rm $curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta
            # Set output variables to input files since no filtering occurred
            output_feature_table=$input_table
            output_fasta=$input_fasta
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

            # For long format: just remove any zero-abundance rows (shouldn't exist, but check)
            checkerror=$(awk '
            BEGIN {
                FS = OFS = "\t"
            }
            NR == 1 {
                # Print header
                print $0
                next
            }
            {
                # For long format: SampleID, OTU, Abundance (3rd column is abundance)
                # Only keep rows with non-zero abundance
                if ($3 > 0) {
                    print $0
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
            output_feature_table="$curated_dir/${input_table_basename%%.txt}_lenFilt.txt"
            output_fasta="$curated_dir/${fasta_basename%%.fasta}_lenFilt.fasta"
        fi
    fi



    # if not tag-jump finterin and no length filtering, use input files as output
    if [[ ( $filter_tag_jumps == "false" || -z $filter_tag_jumps ) && \
          ( -z $min_length_num || $min_length_num == "0" ) && \
          ( -z $max_length_num || $max_length_num == "0" ) ]]; then
        :   
    else
        # Update the arrays with the final output files for this run
        output_feature_tables=("${output_feature_tables[@]/$table_file/$output_feature_table}")
        output_fastas=("${output_fastas[@]/$fasta_file/$output_fasta}")
    fi

done
echo "OTU tables to be merged = ${output_feature_tables[*]}"
echo "Corresponding fasta files = ${output_fastas[*]}"

#######################
# 4. merge OTU tables #
#######################
printf "Merging OTU tables...\n"

# Merge long format tables (simple concatenation and aggregation)
printf "Merging long format OTU tables...\n"
echo "Tables to merge: ${output_feature_tables[*]}"

# Concatenate all tables (first one with header, rest without header)
first_table=true
for table in "${output_feature_tables[@]}"; do
    if $first_table; then
        cat "$table" > "$output_dir/OTU_table_long_merged.temp"
        first_table=false
    else
        tail -n +2 "$table" >> "$output_dir/OTU_table_long_merged.temp"
    fi
done

# Aggregate duplicates by summing abundances (in case same sample+OTU appears in multiple runs)
printf "SampleID\tOTU\tAbundance\n" > "$output_dir/OTU_table_long.txt"
awk '
BEGIN {
    FS = OFS = "\t"
}
NR > 1 {
    key = $1 "\t" $2  # SampleID + OTU
    abundance[key] += $3
}
END {
    for (key in abundance) {
        if (abundance[key] > 0) {
            print key, abundance[key]
        }
    }
}' "$output_dir/OTU_table_long_merged.temp" >> "$output_dir/OTU_table_long.txt"

rm "$output_dir/OTU_table_long_merged.temp"
printf "Merged long format OTU table created: OTU_table_long.txt\n"

# Merge FASTA files
printf "Merging FASTA files...\n"
cat "${output_fastas[@]}" > "$output_dir/merged_OTUs.fasta"
# Remove duplicate sequences from merged FASTA file
printf "Removing duplicate sequences from merged FASTA file...\n"
checkerror=$(seqkit rmdup -w 0 -s "$output_dir/merged_OTUs.fasta" -o "$output_dir/merged_OTUs_dedup.fasta" 2>&1)
check_app_error
mv "$output_dir/merged_OTUs_dedup.fasta" "$output_dir/OTUs.fasta" && rm "$output_dir/merged_OTUs.fasta"

# Update output variables to point to final merged files
output_feature_table="$output_dir/OTU_table_long.txt"
output_fasta="$output_dir/OTUs.fasta"
export output_feature_table
export output_fasta

### Convert OTU table to wide format (optional for large datasets)
# Count unique samples to determine if wide format conversion is feasible
sample_count=$(awk 'NR>1 {samples[$1]=1} END {print length(samples)}' "$output_dir/OTU_table_long.txt")
printf "Number of samples in dataset: $sample_count\n"

# Wide format is problematic with many samples (memory, line length, processing time)
MAX_SAMPLES_FOR_WIDE=5000

if [[ $sample_count -le $MAX_SAMPLES_FOR_WIDE ]]; then
    printf "Converting OTU table to wide format...\n"
    awk '
    BEGIN {
        FS = OFS = "\t"
    }
    NR == 1 {
        # Skip header
        next
    }
    {
        sample = $1
        otu = $2
        abundance = $3
        
        # Store abundance
        data[otu, sample] = abundance
        
        # Track unique OTUs and samples
        if (!(otu in otus)) {
            otus[otu] = ++otu_count
        }
        if (!(sample in samples)) {
            samples[sample] = ++sample_count
            sample_order[sample_count] = sample
        }
    }
    END {
        # Print header: OTU followed by sample names in order they appeared
        printf "OTU"
        for (i = 1; i <= sample_count; i++) {
            printf "\t%s", sample_order[i]
        }
        printf "\n"
        
        # Print data rows
        for (otu in otus) {
            printf "%s", otu
            for (i = 1; i <= sample_count; i++) {
                sample = sample_order[i]
                # Print abundance or 0 if not present
                if ((otu, sample) in data) {
                    printf "\t%s", data[otu, sample]
                } else {
                    printf "\t0"
                }
            }
            printf "\n"
        }
    }' "$output_dir/OTU_table_long.txt" > "$output_dir/OTU_table.txt"
    
    printf "Wide format OTU table created: OTU_table.txt\n"
else
    printf "WARNING: Dataset has $sample_count samples (threshold: $MAX_SAMPLES_FOR_WIDE).\n"
    printf "Skipping wide format conversion:\n"
fi 

### CLEAN UP ###
# Remove all files and folders except output files
printf "Cleaning up temporary files...\n"
if [[ $debugger != "true" ]]; then
    if [[ -f "$output_dir/OTU_table.txt" ]]; then
        # Keep both wide and long format
        find "$output_dir" -mindepth 1 -not \( -name "OTUs.fasta" -o -name "OTU_table_long.txt" -o -name "OTU_table.txt" \) -exec rm -rf {} +
    else
        # Keep only long format (wide format not created)
        find "$output_dir" -mindepth 1 -not \( -name "OTUs.fasta" -o -name "OTU_table_long.txt" \) -exec rm -rf {} +
    fi

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
count_features "$output_dir/OTU_table_long.txt"

end=$(date +%s)
runtime=$((end-start))

printf "# Merged sequencing runs: 

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

 # 1. clustering all samples from all sequencing runs: 
        similarity_threshold = $id
        OTU_type = $otutype
        strands = $strands
        remove_singletons = $remove_singletons
        similarity_type = $simtype
        maxaccepts = $maxaccepts
        mask = $mask
 # 2. Split OTU tables per run\n" > $output_dir/README.txt
 if [[ $curate_otu_table == "true" ]]; then
    printf "# 3. curate OTU tables per run:
    f_value = $f_value
    p_value = $p_value
    $min_length
    $max_length
    
# 4. merge (curated) OTU tables\n" >> $output_dir/README.txt
 else
    printf "# 3. merge OTU tables\n" >> $output_dir/README.txt
 fi

printf "
Output files:
------------
# OTU_table_long.txt = merged OTU abundance table (long format: SampleID, OTU, Abundance)
" >> $output_dir/README.txt

if [[ -f "$output_dir/OTU_table.txt" ]]; then
    printf "# OTU_table.txt      = merged OTU abundance table (wide format: OTU, Sample1, Sample2, ...)\n" >> $output_dir/README.txt
else
    printf "# OTU_table.txt      = NOT CREATED (dataset has $sample_count samples, exceeds $MAX_SAMPLES_FOR_WIDE threshold)\n" >> $output_dir/README.txt
fi

printf "# OTUs.fasta         = merged OTU sequences

Number of OTUs                       = $feature_count
Number of sequences in the OTU table = $nSeqs
Number of samples in the OTU table   = $nSample " >> $output_dir/README.txt

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

#!/bin/bash

# Filter tag-jumps, filter ASV by length, collaplse identical ASVs in DADA2 workflow
# Filter tag-jumps, filter OTUs by length in OTUs workflows
# Input = ASV/OTU table file, txt format, may have "Sequence" column
# 1. tag-jumps filtering, if ON 
# 2. min-max length filtering, if ON
# 3. collapsing identical ASVs, if ON
# 4. only length filtering, if collapseNoMismatch = false

# edit 21.01.2025: multiRunDir handling, 
    # replacing dada2 'collapseNoMismatch' function with usearch_global
    # adding min-max length filtering
##########################################################
###Third-party applications:
#vsearch, R(data.table, ggplot2), awk, seqkit
##########################################################

# Checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
R_version=$(R --version | head -n1 | cut -d " " -f3)
awk_version=$(awk --version | head -n1)
printf "# Checking tool versions ...\n"
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"
printf "# R (version $R_version)\n"
printf "# awk (version $awk_version)\n"

# load variables
f_value=${f_value} # f-value for tag-jumps filtering (UNCROSS2) 
p_value=${p_value} # p-value for tag-jumps filtering (UNCROSS2) 
collapseNoMismatch=${collapseNoMismatch} # collapse identical ASVs (usearch_global --id 1)
max_length=${max_length} # max length of ASVs to keep, if 0 then no max length filter
min_length=${min_length} # min length of ASVs to keep, if 0 then no min length filter

# min_length and max_length handling
if [[ $min_length == "0" ]]; then
    min_length=""
else
    min_length_seqkit="$min_length"
    min_length="--minseqlength $min_length"
fi
if [[ $max_length == "0" ]]; then
    max_length=""
else
    max_length_seqkit="$max_length"
    max_length="--maxseqlength $max_length"
fi

# collapse identical ASVs handling; for OTUs workflow [disabled in OTUs workflow]
if [[ $collapseNoMismatch == "" ]]; then
    collapseNoMismatch="false"
fi

# if f_value == 0 or p_value == 0, then no tag-jumps filtering
if [[ $f_value == "0" ]]; then
    filter_tag_jumps="false"
elif [[ $p_value == "0" ]]; then
    filter_tag_jumps="false"
else
    filter_tag_jumps="true" 
fi

# source for functions
source /scripts/submodules/framework.functions.sh

# check if need to work with multiple or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "CURATE TABLE for multiple sequencing runs in multiRunDir"
    echo "Process = Features (ASVs/OTUs) curation."
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")
    echo "Working in dirs:"
    echo $DIRS
    multiDir="true"
    export multiDir
else
    cd /input
    echo "Working with individual sequencing run"
    echo "Process = Features (ASVs/OTUs) curation."
    DIRS="/input"
    printf "\n workingDirs = $DIRS \n"
    multiDir="false"
    export multiDir
fi

#############################
### Start of the workflow ###
#############################
### looping through multiple sequencing runs (dirs in multiRunDir) 
 # if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)
    cd $seqrun

    # Multi-sequencing run (full pipeline)
    if [[ $multiDir == "true" ]]; then
        # Check for input table; define output_dir and workingDir
        if [[ -f "/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/ASVs_table.txt" ]]; then
            feature_table_file="/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/ASVs_table.txt"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_file="/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/ASVs.fasta"
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/filtered"
            mkdir -p $output_dir
            export output_dir
            workingDir=$output_dir
            printf "\n ASV table filtering, input = $feature_table_file\n"
        elif [[ -f "/input/multiRunDir/${seqrun%%/*}/clustering_out/OTU_table.txt" ]]; then
            feature_table_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/OTU_table.txt"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/OTUs.fasta"
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/multiRunDir/${seqrun%%/*}/clustering_out/filtered"
            mkdir -p $output_dir
            export output_dir
            workingDir=$output_dir
            printf "\n OTU table filtering, input = $feature_table_file\n"
        elif [[ -f "/input/multiRunDir/${seqrun%%/*}/clustering_out/zOTU_table.txt" ]]; then
            feature_table_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/zOTU_table.txt"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/zOTUs.fasta"
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/multiRunDir/${seqrun%%/*}/clustering_out/filtered"
            mkdir -p $output_dir
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n zOTU table filtering, input = $feature_table_file\n"
        else
            printf '%s\n' "ERROR]: Could not find input table.
            Looked for:
            - /ASVs_out.dada2/ASVs_table.txt
            - /clustering_out/OTU_table.txt
            - /clustering_out/zOTU_table.txt
            Please check if ASV/OTU table exists.
            >Quitting" >&2
            end_process
        fi
    # Single sequencing run (full pipeline)
    elif [[ $multiDir == "false" ]]; then   
        # Check for input table; define output_dir and workingDir
        if [[ -f "/input/ASVs_out.dada2/ASVs_table.txt" ]]; then
            feature_table_file="/input/ASVs_out.dada2/ASVs_table.txt"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_file="/input/ASVs_out.dada2/ASVs.fasta"
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/ASVs_out.dada2/filtered"
            mkdir -p $output_dir
            export output_dir
            printf "\n ASV table filtering, input = $feature_table_file\n"
        elif [[ -f "/input/clustering_out/OTU_table.txt" ]]; then
            feature_table_file="/input/clustering_out/OTU_table.txt"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_file="/input/clustering_out/OTUs.fasta"
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/clustering_out/filtered"
            mkdir -p $output_dir
            export output_dir
            printf "\n OTU table filtering, input = $feature_table_file\n"
        elif [[ -f "/input/clustering_out/zOTU_table.txt" ]]; then
            feature_table_file="/input/clustering_out/zOTU_table.txt"
            fasta_file="/input/clustering_out/zOTUs.fasta"
            feature_table_base_name=$(basename $feature_table_file)
            fasta_base_name=$(basename $fasta_file)
            output_dir="/input/clustering_out/filtered"
            mkdir -p $output_dir
            export output_dir
            printf "\n zOTU table filtering, input = $feature_table_file\n"
        else
        printf '%s\n' "ERROR]: CURATE TABLE. Could not find input table.
            Looked for:
            - /ASVs_out.dada2/ASVs_table.txt
            - /clustering_out/OTU_table.txt
            - /clustering_out/zOTU_table.txt
            Please check if ASV/OTU table exists.
            >Quitting" > ERROR_curate_table.log
            printf '%s\n' "ERROR]: CURATE TABLE. Could not find input table.
            Looked for:
            - /ASVs_out.dada2/ASVs_table.txt
            - /clustering_out/OTU_table.txt
            - /clustering_out/zOTU_table.txt
            Please check if ASV/OTU table exists.
            >Quitting" >&2
            end_process
        fi
    fi
 
    # Initialize output variables with input files as defaults at the start
    output_feature_table=$feature_table_file
    output_fasta=$fasta_file

    ### Process samples with UNCROSS2 (tag-jumps filtering) in R
    if [[ $filter_tag_jumps == "true" ]]; then
        printf "# Running tag-jumps filtering (UNCROSS2) for $seqrun\n "
        Rlog=$(Rscript /scripts/submodules/tag_jump_removal.R $feature_table_file $f_value $p_value $fasta_file 2>&1)
        # Check if R script executed successfully
        if [ $? -ne 0 ]; then
            log_error "tag-jumps filtering R script failed with the following error:
            $Rlog
            Please check the parameters and input file.
            >Quitting"
            end_process
        fi
        echo "$Rlog" > "$output_dir/tag-jumps_filt.log"
        wait
    
        # Check if output files were created
        if [ -z "$(find "$output_dir" -name "*_TagJumpFilt.txt")" ]; then
            log_error "UNCROSS2 did not generate the expected output file.
            Please check the log file at $output_dir/tag-jumps_filt.log
            >Quitting"
            end_process
        fi

        printf "\n tag-jumps filtering completed \n"
        # Update output variables if only tag-jumps filtering is performed
        if [[ $collapseNoMismatch != "true" ]] && [[ -z $min_length || $min_length == "0" ]] && [[ -z $max_length || $max_length == "0" ]]; then
            output_feature_table=$output_dir/${feature_table_base_name%%.txt}_TagJumpFilt.txt
            output_fasta=$output_dir/$fasta_base_name
        fi
        # Set input table for next steps
        input_table=$output_dir/${feature_table_base_name%%.txt}_TagJumpFilt.txt
        export input_table
        
        # copy fasta file to output_dir and verify copy succeeded
        if ! cp "$fasta_file" "$output_dir/$fasta_base_name"; then
            log_error "Failed to copy fasta file to output directory
            >Quitting"
            end_process
        fi
        # new fasta_file var in output_dir
        fasta_file=$output_dir/$fasta_base_name
        
        # count ASVs and verify file exists
        if [[ ! -f $fasta_file ]]; then
            log_error "Fasta file not found: $fasta_file
            >Quitting"
            end_process
        fi
        ASVs_count=$(grep -c "^>" "$fasta_file")
        
        # format R-log file
        sed -i "s/;; /\n/g" "$output_dir/tag-jumps_filt.log"

    elif [[ $filter_tag_jumps == "false" ]]; then
        printf "# Skipping tag-jumps filtering\n"
        input_table=$feature_table_file
        export input_table
        # fasta_file var does not change here
    fi

    ###################################
    ### Filtering the ASV/OTU table ###
    ###################################
    # collapse identical ASVs and length filtering
    if [[ $collapseNoMismatch == "true" ]]; then
        printf "Starting collapseNoMismatch (with length filter if ON) ... \n"
        # count input ASVs
        ASVs_count=$(grep -c "^>" $fasta_file)

        # dereplicate
        checkerror=$(vsearch --derep_fulllength $fasta_file \
                --output ${fasta_file%%.fasta}_derep.fasta \
                --fasta_width 0 \
                --uc ${fasta_file%%.fasta}_derep.uc \
                $min_length $max_length 2>&1)
        check_app_error

        # count ASVs after dereplication
        ASVs_derep_count=$(grep -c "^>" ${fasta_file%%.fasta}_derep.fasta)

        # skip collapsing and length filtering if ASVs_derep_count is 0
        if [[ $ASVs_derep_count == "0" ]]; then
            ASVs_collapsed_result="All Features (ASVs/OTUs) were filtered out based on the length filter 
            (min_length $min_length bp and max_length $max_length bp).
            Skipping collapseNoMismatch and length filtering. "
            echo $ASVs_collapsed_result
            # Set output variables to input files since no filtering occurred
            output_feature_table=$feature_table_file
            output_fasta=$fasta_file
            continue
        else
            ASVs_lenFilt_result="ASVs after length filtering = $ASVs_derep_count (input was $ASVs_count ASVs)"
            echo $ASVs_lenFilt_result
        fi

        # usearch_global to find exact matches (no mismatch), excluding terminal gaps
        checkerror=$(vsearch --usearch_global ${fasta_file%%.fasta}_derep.fasta \
                --db ${fasta_file%%.fasta}_derep.fasta \
                --id 1 \
                --uc $output_dir/usearch_global.uc \
                --strand plus \
                --maxaccepts 99999 2>&1)
        check_app_error
       
        ### Convert wide format OTU table to long format in awk (bash) 
        #  handles "Sequence" column
        checkerror=$(cat $input_table \
        | awk '
        BEGIN { 
        FS="\t"; OFS="\t";
        print "OTU", "SampleID", "Abundance"   # Header of the resulting table
        }
        NR==1 {
        seq_col = -1;                      # Initialize sequence column index
        for (i=2; i<=NF; i++) {
            if ($i == "Sequence") {        # Check if column is "Sequence"
            seq_col = i;                   # Store sequence column index
            continue;                      # Skip this column
            }
            sampleIDs[i] = $i;             # Store sample IDs from the header row
        }
        }
        NR>1 {
        otu = $1;                          # Get the OTU ID from the first column
        for (i=2; i<=NF; i++) {
            if (i == seq_col) continue;    # Skip the sequence column
            if ($i > 0) {                  # Skip zero abundances
            print otu, sampleIDs[i], $i;   # Print OTU ID, SampleID, and Abundance
            }
        }
        }' > $output_dir/table_long.txt 2>&1)
        check_app_error

        # make collapsed ASV table
        Rlog=$(Rscript /scripts/submodules/make_FeatureTable_collapseASVs.R \
                --derepuc ${fasta_file%%.fasta}_derep.uc \
                --uc $output_dir/usearch_global.uc \
                --asv $output_dir/table_long.txt \
                --output $output_dir/${feature_table_base_name%%.txt}_collapsed.txt)
        # Check if R script executed successfully
        if [ $? -ne 0 ]; then
            log_error "make_FeatureTable_collapseASVs.R script failed with the following error:
            $Rlog
            >Quitting"
            end_process
        fi
        echo $Rlog > $output_dir/make_FeatureTable_collapseASVs.log
        wait    
        # format R-log file
        sed -i "s/;; /\n/g" $output_dir/make_FeatureTable_collapseASVs.log 

        ### extract representative sequences for collapsed ASVs
        # extract first column, skip header
        cut -f1 $output_dir/${feature_table_base_name%%.txt}_collapsed.txt | tail -n +2 > $output_dir/seq_IDs.txt 

        # extract sequences from ASVs_derep.fasta that match the OTU_ids.txt
        seqkit grep -f $output_dir/seq_IDs.txt $fasta_file \
                        -w 0 > $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta

        # count collapsed (and length filtered) ASVs
        if [[ -f $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta ]]; then
            ASVs_collapsed_count=$(grep -c "^>" $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta)
            # for the report, if ASVs_collapsed_count == ASVs_count, then no new files outputted
            if [[ $ASVs_collapsed_count == $ASVs_count ]]; then
                ASVs_collapsed_result=$"Output has the same number of Features (ASVs/OTUs) as input.
Input had $ASVs_count Features; lenFilt = $ASVs_derep_count Features).
No new files outputted."
                echo -e "$ASVs_collapsed_result"
                # remove intermediate files
                rm $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta
                rm $output_dir/${feature_table_base_name%%.txt}_collapsed.txt
                
                # Set output variables based on whether tag-jumps filtering was performed
                if [[ $filter_tag_jumps == "true" ]]; then
                    output_feature_table=$output_dir/${feature_table_base_name%%.txt}_TagJumpFilt.txt
                    output_fasta=$output_dir/$fasta_base_name
                else
                    output_feature_table=$feature_table_file
                    output_fasta=$fasta_file
                fi
            fi
            # for the report, ASVs_collapsed_count < $ASVs_count
            if [[ $ASVs_collapsed_count < $ASVs_count ]]; then
                # output files variables for the Merge sequencing runs
                output_feature_table=$output_dir/${feature_table_base_name%%.txt}_collapsed.txt
                output_fasta=$output_dir/${fasta_base_name%%.fasta}_collapsed.fasta
                # for the report
                ASVs_collapsed_result="Outputted $ASVs_collapsed_count Features (ASVs/OTUs) (lenFilt resulted in $ASVs_derep_count Features).
Input had $ASVs_count Features.
    - ${feature_table_base_name%%.txt}_collapsed.txt: Feature table after collapsing identical Features (ASVs/OTUs). Contains $ASVs_collapsed_count Features.
    - ${fasta_base_name%%.fasta}_collapsed.fasta: Representative sequences after collapsing"
                echo -e "$ASVs_collapsed_result"
            fi
        fi 

    #############################
    ### only length filtering ###
    #############################
    elif [[ $collapseNoMismatch == "false" ]] && \
        [[ $min_length != "0" && -n $min_length ]] && \
        [[ $max_length != "0" && -n $max_length ]]; then
        printf "Filtering by length, min_length = $min_length_seqkit, max_length = $max_length_seqkit. \n"
        
        # get basenames for correct naming of output files
        input_table_base_name=$(basename $input_table)
        fasta_base_name=$(basename $fasta_file)
        # count input ASVs
        ASVs_count=$(grep -c "^>" $fasta_file)
        echo "Feature (ASVs/OTUs) count = $ASVs_count"
        # filter by length
        checkerror=$(seqkit seq -w 0 -g \
                    -m $min_length_seqkit \
                    -M $max_length_seqkit \
                    $fasta_file \
                    > $output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta 2>&1)
        check_app_error

        # count length filtered ASVs and proceed with the rest of the steps
        ASVs_lenFilt=$(grep -c "^>" $output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta)
        echo "length filtered Feature (ASV/OTU) count = $ASVs_lenFilt"
        if [[ $ASVs_lenFilt == 0 ]]; then
            ASVs_lenFilt_result=$"All Features (ASVs/OTUs) were filtered out based on the length filter
            (min_length $min_length_seqkit bp and max_length $max_length_seqkit bp).
            No new files generated.
            Input table was $input_table_base_name and input fasta was $fasta_base_name with $ASVs_count sequences"
            echo -e "$ASVs_lenFilt_result"
            rm $output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta
            # Set output variables to input files since no filtering occurred
            output_feature_table=$input_table
            output_fasta=$fasta_file
            
        elif [[ $ASVs_lenFilt == $ASVs_count ]]; then
            ASVs_lenFilt_result=$"None of the Features (ASVs/OTUs) were filtered out based on the length filter
            (min_length $min_length_seqkit bp and max_length $max_length_seqkit bp).
            No new files generated.
            Input table was $input_table_base_name and input fasta was $fasta_base_name with $ASVs_count sequences"
            echo -e "$ASVs_lenFilt_result"
            export ASVs_lenFilt_result
            rm $output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta
            # Set output variables to input files since no filtering occurred
            output_feature_table=$input_table
            output_fasta=$fasta_file
        else          
            # filter the table
            checkerror=$(seqkit seq \
                        -n $output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta \
                        > $output_dir/${fasta_base_name%%.fasta}_IDs.txt 2>&1)
            check_app_error
            checkerror=$(grep -f $output_dir/${fasta_base_name%%.fasta}_IDs.txt $input_table \
                                > $output_dir/${input_table_base_name%%.txt}_lenFilt.temp 2>&1)
            check_app_error
            # remove intermediate files
            rm $output_dir/${fasta_base_name%%.fasta}_IDs.txt
            # add 1st row of the $input_table to the $output_dir/${input_table_base_name%%.txt}_lenFilt.temp
            header=$(head -n 1 $input_table)
            sed -i "1i\\$header" "$output_dir/${input_table_base_name%%.txt}_lenFilt.temp"

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
            ' "$output_dir/${input_table_base_name%%.txt}_lenFilt.temp" \
            > "$output_dir/${input_table_base_name%%.txt}_lenFilt.txt" 2>&1)
            check_app_error
            # remove intermediate files
            if [[ -f $output_dir/${input_table_base_name%%.txt}_lenFilt.temp ]]; then
                rm $output_dir/${input_table_base_name%%.txt}_lenFilt.temp
            fi
            # output files variables for the Merge sequencing runs
            output_feature_table=$output_dir/${input_table_base_name%%.txt}_lenFilt.txt
            output_fasta=$output_dir/${fasta_base_name%%.fasta}_lenFilt.fasta
            # for the report
            ASVs_lenFilt_result=$"Features (ASVs/OTUs) after length filtering = $ASVs_lenFilt.
Input had $ASVs_count Features.
    - ${input_table_base_name%%.txt}_lenFilt.txt: Feature table after length filtering. Contains $ASVs_lenFilt Features.
    - ${fasta_base_name%%.fasta}_lenFilt.fasta: Representative sequences file after length filtering"
            echo -e "$ASVs_lenFilt_result"
        fi
    fi

    # Remove tempdir2 if debugger is not true
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        if [[ -f $output_dir/make_FeatureTable_collapseASVs.log ]]; then
            rm $output_dir/make_FeatureTable_collapseASVs.log
        fi
        if [[ -f $output_dir/seq_IDs.txt ]]; then
            rm $output_dir/seq_IDs.txt
        fi
        if [[ -f $output_dir/table_long.txt ]]; then
            rm $output_dir/table_long.txt
        fi
        if [[ -f $output_dir/${fasta_base_name%%.fasta}_derep.uc ]]; then
            rm $output_dir/${fasta_base_name%%.fasta}_derep.uc
            rm $output_dir/${fasta_base_name%%.fasta}_derep.fasta
        fi
        if [[ -f $output_dir/usearch_global.uc ]]; then
            rm $output_dir/usearch_global.uc
        fi
    fi

    # Make README.txt file
    end=$(date +%s)
    runtime=$((end-start))
    readme_table_filtering $output_dir $runtime
    
    ### if working with multiRunDir then cd /input/multiRunDir to start next seqrun
        # Store outputs for this sequencing run if working with multiRunDir
    if [[ $multiDir == "true" ]]; then
        # Create arrays if they don't exist yet
        declare -a output_feature_tables
        declare -a output_fastas
        # Add current run's outputs to arrays
        output_feature_tables+=("$output_feature_table")
        output_fastas+=("$output_fasta")
        cd /input/multiRunDir
    else
        cd /input
    fi
done

# output_feature_tables to output_feature_table when multiDir is true
if [[ $multiDir == "true" ]]; then
    output_feature_table=$(IFS=,; echo "${output_feature_tables[*]}")
    output_fasta=$(IFS=,; echo "${output_fastas[*]}")
else
    echo "output_feature_table=$output_feature_table"
    echo "output_fasta=$output_fasta"
fi

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"

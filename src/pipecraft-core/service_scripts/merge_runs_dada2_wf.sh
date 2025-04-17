#!/bin/bash

# Merge sequencing runs processed with DADA2 if working with multuple runs in multiRunDir. 
 # Samples with the same name across runs are merged together.

################################################
###Third-party applications:
# dada2, R
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# DADA2 version: $dada2_version\n"
printf "# vsearch version: $vsearch_version\n"

start_time=$(date)
start=$(date +%s)
# source for functions
source /scripts/submodules/framework.functions.sh

# Excecute only if multiDir = true
if [[ ! -d "/input/multiRunDir" ]]; then
    printf '%s\n' "ERROR]: multiRunDir not detected. Cannot merge sequencing runs. 
    >DONE." >&2
    end_process
elif [[ $merge_runs == "true" ]]; then
    printf "Starting merge sequencing runs...\n"
    #output dir
    output_dir=$"/input/multiRunDir/merged_runs"
    export output_dir
    # remove output dir if it already exists
    if [[ -d "$output_dir" ]]; then
        rm -rf $output_dir
    fi
    # create new output dir
    mkdir -p $output_dir

    echo "input tables: $output_feature_table" 
    echo "input fasta: $output_fasta"
    echo "output dir: $output_dir"
else
    printf '%s\n' "ERROR]: Merge sequencing runs is not enabled. Exiting.\n" >&2
    end_process
fi

# Checking tool versions
R_version=$(R --version | head -n1 | cut -d " " -f3)
dada2_version=$(Rscript -e "packageVersion('dada2')" 2>/dev/null | awk '{print $2}' | sed -e "s/‘//g" -e 's/’//g')
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# Checking tool versions ...\n"
printf "# R (version $R_version)\n"
printf "# DADA2 (version $dada2_version)\n"
printf "# vsearch (version $vsearch_version)\n"

# load variables
collapseNoMismatch=${collapseNoMismatch} # collapse identical ASVs (usearch_global --id 1)

### Merge ASV tables with dada2 mergeSequenceTables function in R
printf "# Running DADA2 mergeSequenceTables ...\n"
Rlog=$(Rscript /scripts/submodules/dada2_mergeRuns.R 2>&1)
echo $Rlog > $output_dir/dada2_mergeRuns.log 
wait
# format R log
sed -i 's/;; /\n/g' $output_dir/dada2_mergeRuns.log

# check if output files exist
if [[ ! -f $output_dir/ASVs.fasta ]] || [[ ! -f $output_dir/ASVs_table.txt ]]; then
    printf '%s\n' "ERROR]: Output files not found. Merge runs FAILED. 
    >Exiting." >&2
    end_process
else
    printf "Merge runs completed successfully. \n"
    fasta_file=$output_dir/ASVs.fasta
    table_file=$output_dir/ASVs_table.txt
    fasta_base_name=${fasta_file##*/}
    feature_table_base_name=${table_file##*/}
fi

###############################################
### collapse identical Features (ASVs/OTUs) ###
###############################################
if [[ $collapseNoMismatch == "true" ]]; then
    printf "Starting collapseNoMismatch ... \n"
    # count input ASVs
    ASVs_count=$(grep -c "^>" $fasta_file)

    # dereplicate
    checkerror=$(vsearch --derep_fulllength $fasta_file \
            --output ${fasta_file%%.fasta}_derep.fasta \
            --fasta_width 0 \
            --uc ${fasta_file%%.fasta}_derep.uc 2>&1)
    check_app_error

    # usearch_global to find exact matches (no mismatch), excluding terminal gaps
    checkerror=$(vsearch --usearch_global ${fasta_file%%.fasta}_derep.fasta \
            --db ${fasta_file%%.fasta}_derep.fasta \
            --id 1 \
            --uc $output_dir/usearch_global.uc \
            --strand plus \
            --maxaccepts 999999 2>&1)
    check_app_error
    
    ### Convert wide format OTU table to long format in awk (bash) 
    #  handles "Sequence" column
    checkerror=$(cat $table_file \
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
    echo "make collapsed ASV table"
    Rlog=$(Rscript /scripts/submodules/make_FeatureTable_collapseASVs.R \
            --derepuc ${fasta_file%%.fasta}_derep.uc \
            --uc $output_dir/usearch_global.uc \
            --asv $output_dir/table_long.txt \
            --fasta $fasta_file \
            --output $output_dir/${feature_table_base_name%%.txt}_collapsed.txt)
    # Check if R script executed successfully
    if [ $? -ne 0 ]; then
        log_error "make_FeatureTable_collapseASVs.R script failed with the following error:
        $Rlog
        >Quitting"
        end_process
    else 
        rm $output_dir/table_long.txt
        rm $output_dir/usearch_global.uc
        rm ${fasta_file%%.fasta}_derep.uc
        rm ${fasta_file%%.fasta}_derep.fasta
    fi
    echo $Rlog > $output_dir/make_FeatureTable_collapseASVs.log
    wait    
    # format R-log file
    sed -i "s/;; /\n/g" $output_dir/make_FeatureTable_collapseASVs.log 

    ### extract representative sequences for collapsed ASVs
    # extract first column, skip header
    cut -f1 $output_dir/${feature_table_base_name%%.txt}_collapsed.txt | tail -n +2 > $output_dir/seq_IDs.txt 

    # extract sequences from ASVs_derep.fasta that match the OTU_ids.txt
    checkerror=$(seqkit grep -f $output_dir/seq_IDs.txt $fasta_file \
                    -w 0 > $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta 2>&1)
    check_app_error
    rm $output_dir/seq_IDs.txt
    # count collapsed (and length filtered) ASVs
    if [[ -f $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta ]]; then
        ASVs_collapsed_count=$(grep -c "^>" $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta)
        # for the report, if ASVs_collapsed_count == ASVs_count, then no new files outputted
        if [[ $ASVs_collapsed_count == $ASVs_count ]]; then
            ASVs_collapsed_result=$"Output has the same number of Features (ASVs/OTUs) as input. No new files outputted."
            echo -e "$ASVs_collapsed_result"
            # remove intermediate files
            rm $output_dir/${fasta_base_name%%.fasta}_collapsed.fasta
            rm $output_dir/${feature_table_base_name%%.txt}_collapsed.txt
           
        fi
        # for the report, ASVs_collapsed_count < $ASVs_count
        if [[ $ASVs_collapsed_count < $ASVs_count ]]; then
            # for the report
            ASVs_collapsed_result="Outputted $ASVs_collapsed_count Features (ASVs/OTUs).
Input had $ASVs_count Features.
- ${feature_table_base_name%%.txt}_collapsed.txt = Feature table after collapsing identical Features (ASVs/OTUs). Contains $ASVs_collapsed_count Features.
- ${fasta_base_name%%.fasta}_collapsed.fasta = Representative sequences after collapsing"
            echo -e "$ASVs_collapsed_result"
        fi
    fi 
fi

# Remove intermediate files if debugger is not true
if [[ $debugger != "true" ]]; then
    rm $output_dir/make_FeatureTable_collapseASVs.log 
    rm $output_dir/dada2_mergeRuns.log
    if [[ -f "/input/multiRunDir/.curate_table_params" ]]; then
        rm /input/multiRunDir/.curate_table_params
    fi
fi

# count features and sequences; outputs variables feature_count, nSeqs, nSample
count_features "$output_dir/ASVs_table.txt"

### Make README.txt file (merged_runs)
end=$(date +%s)
runtime=$((end-start))
printf "# Merged sequencing runs with DADA2 mergeSequenceTables function.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input tables:\n" > $output_dir/README.txt

# Add each input table path
IFS=',' read -ra TABLES <<< "$output_feature_table"
for table in "${TABLES[@]}"; do
    printf "%s\n" "$table" >> $output_dir/README.txt
done

printf "
Output files:
------------
# ASVs_table.txt = merged ASV abundance table
# ASVs.fasta     = merged ASV sequences

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample " >> $output_dir/README.txt


if [[ $collapseNoMismatch == "true" ]]; then
    printf "\n
Outputs after CollapsedNoMismatch = true:
-----------------------------------------
$ASVs_collapsed_result \n" >> $output_dir/README.txt
fi

if [[ $collapseNoMismatch == "true" ]] && [[ -f $output_dir/${feature_table_base_name%%.txt}_collapsed.txt ]]; then
    count_features $output_dir/${feature_table_base_name%%.txt}_collapsed.txt
    printf "
${feature_table_base_name%%.txt}_collapsed.txt = merged and collapsed ASV abundance table
${fasta_base_name%%.fasta}_collapsed.fasta = merged and collapsed ASV sequences

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample " >> $output_dir/README.txt
fi

printf "\n
#############################################
###Third-party applications for this process:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
# vsearch (version $vsearch_version)
    #citation: Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584
    #https://github.com/torognes/vsearch
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

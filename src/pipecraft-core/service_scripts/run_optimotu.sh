#!/bin/bash

# OptimOTU workflow


# Start time
start_time=$(date)
start=$(date +%s)

# Copy the configuration file
cp -f /scripts/pipeline_options.yaml /optimotu_targets/pipeline_options.yaml # to targets
cp -f /scripts/pipeline_options.yaml /optimotu_targets/sequences/pipeline_options.yaml # to user

# Activate the conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate OptimOTU_v4

# Print current environment for debugging
echo "PATH: $PATH"
echo "Current directory: $(pwd)"
echo "Conda environment: $CONDA_DEFAULT_ENV"

# Run the OptimOTU pipeline
Rlog=$(R -e "targets::tar_config_set(script='/optimotu_targets/_targets.R'); targets::tar_make(callr_function=NULL)" 2>&1)
echo $Rlog >> /optimotu_targets/sequences/optimotu_targets.log
wait

# Check if the R command was successful
if [ $? -eq 0 ]; then
    echo "OptimOTU pipeline completed successfully."
    echo "OptimOTU pipeline completed successfully." >> /optimotu_targets/sequences/optimotu_targets.log

    # Copy the output folder to the sequences folder
    echo "Copying output files to sequences folder..."
    cp -r /optimotu_targets/output/* /optimotu_targets/sequences/
    
    # Check if the copy was successful
    if [ $? -eq 0 ]; then
        echo "Files successfully copied to /optimotu_targets/sequences/"
        echo "Files successfully copied to /optimotu_targets/sequences/" >> /optimotu_targets/sequences/optimotu_targets.log
    else
        echo "Error: Failed to copy files to /optimotu_targets/sequences/"
        echo "Error: Failed to copy files to /optimotu_targets/sequences/" >> /optimotu_targets/sequences/optimotu_targets.log
        echo "Copy error details:" >> /optimotu_targets/sequences/optimotu_targets.log
        cp -r /optimotu_targets/output/* /optimotu_targets/sequences/ 2>> /optimotu_targets/sequences/optimotu_targets.log
        exit 1
    fi
else
    echo "Error: OptimOTU pipeline failed."
    echo "Error: OptimOTU pipeline failed." >> /optimotu_targets/sequences/optimotu_targets.log
    exit 1
fi

echo "All operations completed."

if [ ! -z "$HOST_UID" ] && [ ! -z "$HOST_GID" ]; then
  echo "Setting ownership of /optimotu_targets to $HOST_UID:$HOST_GID"
  chown -R $HOST_UID:$HOST_GID /optimotu_targets/sequences
fi

### Make README.txt 
end=$(date +%s)
runtime=$((end-start))
in=$(echo $in_name | sed -e 's/\/extraFiles\///')
printf "# OptimOTU workflow:

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

The outputs of the pipeline are a set of tables in TSV format (tab-delimited files) and 
RDS format (for easy loading in R), as well as sequences in gzipped FASTA format.

The plausible and reliable ({conf}) versions of the table are based on taxonomic assignments at the 50% and 90% probability thresholds, respectively.

Output files:
# asv_table               = ASV table as a sparse matrix (long format) with five columns: sample, seqrun, seq_id, seq_idx, and nread.
# asv2tax_{conf}          = Taxonomic assignments for each ASV at the 50% and 90% probability thresholds, respectively.
# otu_taxonomy_{conf}     = Taxonomy for each OTU at the 50% and 90% probability thresholds, respectively.
# otu_table_sparse_{conf} = OTU table as a sparse matrix (long format) with five columns: sample, seqrun, seq_id, seq_idx, and nread.
# otu_table_{conf}        = OTU table as a dense matrix (wide format) with columns as samples and rows as OTUs.
# otu_{conf}.fasta        = representative OTU sequences for the 50% and 90% probability thresholds for plausible and reliable OTUs, respectively.
# read_counts_{conf}.tsv  = the number of reads in each sample present after each stage of the pipeline.
# optimotu_targets.log    = R log file about the OptimOTU pipeline

All output files are also zipped into OptimOTU_in_PipeCraft2_*.zip (except for the log file).

##############################################
###Third-party applications for this process:
# OptimOTU pipeline v5.0.0 (https://github.com/brendanf/optimotu_targets/releases/tag/v5.0.0)
  # wraps optimotu=0.9.3 and optimotu.pipeline=0.5.2
    citation: Furneaux, B., Anslan, S., Somervuo, P., Hultman, J., Abrego, N., Roslin, T., & Ovaskainen, O. (2025). OptimOTU: Taxonomically aware OTU clustering with optimized thresholds and a bioinformatics workflow for metabarcoding data. arXiv preprint arXiv:2502.10350.
##############################################" > /optimotu_targets/sequences/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "
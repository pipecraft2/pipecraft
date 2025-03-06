#!/bin/bash

# Copy the configuration file
cp -f /scripts/pipeline_options.yaml /optimotu_targets/pipeline_options.yaml

# Activate the conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate OptimOTU_v4

# Print current environment for debugging
echo "PATH: $PATH"
echo "Current directory: $(pwd)"
echo "Conda environment: $CONDA_DEFAULT_ENV"

# Run the OptimOTU pipeline
R -e "targets::tar_config_set(script='/optimotu_targets/_targets.R'); targets::tar_make(callr_function=NULL)"

# Check if the R command was successful
if [ $? -eq 0 ]; then
    echo "OptimOTU pipeline completed successfully."

    
    # Copy the output folder to the sequences folder
    echo "Copying output files to sequences folder..."
    cp -r /optimotu_targets/output/* /optimotu_targets/sequences/
    
    # Check if the copy was successful
    if [ $? -eq 0 ]; then
        echo "Files successfully copied to /optimotu_targets/sequences/"
    else
        echo "Error: Failed to copy files to /optimotu_targets/sequences/"
        exit 1
    fi
else
    echo "Error: OptimOTU pipeline failed."
    exit 1
fi

echo "All operations completed."

if [ ! -z "$HOST_UID" ] && [ ! -z "$HOST_GID" ]; then
  echo "Setting ownership of /optimotu_targets to $HOST_UID:$HOST_GID"
  chown -R $HOST_UID:$HOST_GID /optimotu_targets/sequences
fi
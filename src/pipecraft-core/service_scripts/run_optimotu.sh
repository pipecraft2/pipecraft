#!/bin/bash

# Set proper encoding for output
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export LC_CTYPE=C.UTF-8

#####################
# OptimOTU workflow #
#####################

# # Check if input files match expected format
fileFormat=${fileFormat}
rawFilesDir=${rawFilesDir}
export fileFormat

echo "specified fileFormat: $fileFormat"
echo "specified rawFilesDir: $rawFilesDir"

# Verify 01_raw exists (from bind). Do not create or modify it here.
RAW_BASE="/optimotu_targets/sequences"
LINK_DIR="${RAW_BASE}/01_raw"

if [ ! -d "${LINK_DIR}" ]; then
  echo "[ERROR]: Required directory not found: ${LINK_DIR}. Ensure the selected runs directory is bind-mounted to /optimotu_targets/sequences/01_raw" >&2
  echo "[ERROR]: Missing required directory: ${LINK_DIR}" > "${RAW_BASE}/optimotu_targets.log"
  exit 1
fi

ls -la "${RAW_BASE}"
# # check files in 01_raw and subdirectories
# # Find all subdirectories in 01_raw
 while IFS= read -r subdir; do
     echo "Checking files in $subdir..."
    
     # Count files with the specified extension in the subdirectory
     # Follow symlinks when counting
     file_count=$(find -L "$subdir" -maxdepth 1 -type f -name "*.${fileFormat}" | wc -l)
    
     if [ "$file_count" -eq 0 ]; then
         printf "[ERROR]: No %s files found in %s\n" "${fileFormat}" "${subdir}" >&2
         printf "[ERROR]: No %s files found in %s\n" "${fileFormat}" "${subdir}" > /optimotu_targets/sequences/optimotu_targets.log
         exit 1
     else
         echo "Found $file_count .${fileFormat} files in $subdir"
     fi
done < <(find -L /optimotu_targets/sequences/01_raw -mindepth 1 -maxdepth 1 -type d ! -name 01_raw)

echo "All directories contain valid .${fileFormat} files. Proceeding with pipeline."

# Start time
start_time=$(date)
start=$(date +%s)
echo $fileFormat
echo $readType

# Copy the configuration file
cp -f /scripts/pipeline_options.yaml /optimotu_targets/pipeline_options.yaml # to targets
cp -f /scripts/pipeline_options.yaml /optimotu_targets/sequences/pipeline_options.yaml # to user

# Activate the conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate OptimOTU_v5

# Try to install a compatible version of qs2 only on ARM64 macOS
echo "Checking system architecture and OS..."
echo "Host OS: ${HOST_OS:-unknown}"
echo "Host Architecture: ${HOST_ARCH:-unknown}"

if [ "${HOST_OS}" = "mac" ] && [ "${HOST_ARCH}" = "arm64" ]; then
    echo "ARM64 macOS detected. Installing compatible version of qs2..."
    R --vanilla -e '
      # Remove existing qs2 if present
      if("qs2" %in% installed.packages()[,"Package"]) remove.packages("qs2")
      
      # Install from source with minimal optimizations
      Sys.setenv(PKG_CFLAGS="-O0 -march=x86-64")
      Sys.setenv(PKG_CXXFLAGS="-O0 -march=x86-64")
      install.packages("qs2", type="source", repos="https://cloud.r-project.org")
    '
else
    echo "Not ARM64 macOS. Skipping qs2 installation."
fi


cd /optimotu_targets

# Print current environment for debugging
echo "PATH: $PATH"
echo "Current directory: $(pwd)"
echo "Conda environment: $CONDA_DEFAULT_ENV"

# Run the OptimOTU pipeline and capture output with live updates
echo "Starting OptimOTU pipeline..."
echo "Starting OptimOTU pipeline..." >> /optimotu_targets/sequences/optimotu_targets.log

# Use tee to show live output while capturing to log
R --vanilla -e "targets::tar_config_set(script='/optimotu_targets/_targets.R'); targets::tar_make(callr_function=NULL)" 2>&1 | tee -a /optimotu_targets/sequences/optimotu_targets.log
R_EXIT_STATUS=${PIPESTATUS[0]}


# Check the actual R script exit status
if [ $R_EXIT_STATUS -eq 0 ]; then
    echo "OptimOTU pipeline completed successfully."
else
    echo "Error: OptimOTU pipeline failed with exit status $R_EXIT_STATUS"
    echo "Error: OptimOTU pipeline failed with exit status $R_EXIT_STATUS" >> /optimotu_targets/sequences/optimotu_targets.log
    echo "Last 50 lines of R output from log file:"
    tail -n 50 /optimotu_targets/sequences/optimotu_targets.log
fi

# Copy the output folder to the sequences folder, regardless of success/failure
echo "Copying output files to sequences folder..."
if cp -r /optimotu_targets/output/* /optimotu_targets/sequences 2>/dev/null; then
    echo "Files successfully copied to /optimotu_targets/sequences/"
    echo "Files successfully copied to /optimotu_targets/sequences/" >> /optimotu_targets/sequences/optimotu_targets.log
else
    echo "Warning: Failed to copy files to /optimotu_targets/sequences/"
    echo "Warning: Failed to copy files to /optimotu_targets/sequences/" >> /optimotu_targets/sequences/optimotu_targets.log
    echo "Copy error details:" >> /optimotu_targets/sequences/optimotu_targets.log
    cp -r /optimotu_targets/output/* /optimotu_targets/sequences/ 2>> /optimotu_targets/sequences/optimotu_targets.log 2>/dev/null || echo "Fallback copy also failed"
fi

# Exit with R status after copying output
if [ $R_EXIT_STATUS -ne 0 ]; then
    exit $R_EXIT_STATUS
fi


### Make README.txt 
end=$(date +%s)
runtime=$((end-start))

cat << EOF > /optimotu_targets/sequences/01_raw/README.txt
# OptimOTU workflow:

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

The outputs of the pipeline are a set of tables in TSV format (tab-delimited files) and 
RDS format (for easy loading in R), as well as sequences in gzipped FASTA format.

The plausible and reliable (conf) versions of the table are based on taxonomic assignments at the 50% and 90% probability thresholds, respectively.

Output files:
# asv_table               = ASV table as a sparse matrix (long format) with five columns: sample, seqrun, seq_id, seq_idx, and nread.
# asv2tax_(conf)          = Taxonomic assignments for each ASV at the 50% and 90% probability thresholds, respectively.
# otu_taxonomy_(conf)     = Taxonomy for each OTU at the 50% and 90% probability thresholds, respectively.
# otu_table_sparse_(conf) = OTU table as a sparse matrix (long format) with five columns: sample, seqrun, seq_id, seq_idx, and nread.
# otu_table_(conf)        = OTU table as a dense matrix (wide format) with columns as samples and rows as OTUs.
# otu_(conf).fasta        = representative OTU sequences for the 50% and 90% probability thresholds for plausible and reliable OTUs, respectively.
# read_counts_(conf).tsv  = the number of reads in each sample present after each stage of the pipeline.
# optimotu_targets.log    = R log file about the OptimOTU pipeline

All output files are also zipped into OptimOTU_in_PipeCraft2_*.zip (except for the log file).

##############################################
###Third-party applications for this process:
# OptimOTU pipeline v5.0.0 (https://github.com/brendanf/optimotu_targets/releases/tag/v5.0.0)
    citation: Furneaux, B., Anslan, S., Somervuo, P., Hultman, J., Abrego, N., Roslin, T., & Ovaskainen, O. (2025). OptimOTU: Taxonomically aware OTU clustering with optimized thresholds and a bioinformatics workflow for metabarcoding data. arXiv preprint arXiv:2502.10350.
##############################################
EOF

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "
echo "All operations completed."


if [ ! -z "$HOST_UID" ] && [ ! -z "$HOST_GID" ]; then
  echo "Setting ownership of /optimotu_targets to $HOST_UID:$HOST_GID"
  chown -R $HOST_UID:$HOST_GID /optimotu_targets/sequences
fi
#!/bin/bash
  # $1 = name of the directory containing FASTQ files
  echo -e "\n"
  echo -e "Input data: " $1
  echo -e "Output: " Step1_Results/"$1"
  echo -e "Temporary workdirs: " Step1_WorkDirs/"$1"
  BASEDIR=$(pwd)
  ## Create output directories
  mkdir -p Step1_Results/"$1"
  mkdir -p Step1_WorkDirs/"$1"
  cd Step1_WorkDirs/"$1"


## Step-1 - with pre-demultiplexed data
stdbuf -oL -eL nextflow run \
  vmikk/NextITS -r main \
  -resume \
  --step "Step1" \
  --storagemode "copy" \
  -params-file /scripts/NextFlowConfig.json \
  --input      "$BASEDIR"Input/"$1" \
  --outdir     "$BASEDIR"Input/Step1_Results/"$1" \
  -work-dir    "$BASEDIR"Input/Step1_WorkDirs/"$1" \
  --tracedir   "$BASEDIR"Input/Step1_WorkDirs/"$1"/pipeline_info \
  --demultiplexed true \
  -ansi-log    false \
  2>&1 | tee -a "$BASEDIR"/Input/Nextflow__"$1".log
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
  ## Run Step-1 of the pipeline
  nextflow run /scripts/NextITS/main.nf \
    -resume \
    -params-file /scripts/NextFlowConfig.json \
    --demultiplexed true \
    --input      "$BASEDIR"/Input/"$1" \
    --outdir     "$BASEDIR"/Step1_Results/"$1" \
    -work-dir    "$BASEDIR"/Step1_WorkDirs/"$1" \
    --tracedir   "$BASEDIR"/Step1_WorkDirs/"$1"/pipeline_info \
    -ansi-log    false \
    >> Nextflow__"$1".log

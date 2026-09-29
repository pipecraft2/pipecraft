#!/bin/bash
  # $1 = name of the directory containing FASTQ files
  export NXF_HOME="${NXF_HOME:-/Input/.nextflow}"
  mkdir -p "$NXF_HOME"
  export NXF_ANSI_LOG="false"
  export NXF_LOG_COLOR="false"
  export NXF_ANSI="false"
  export TERM="dumb"
  echo -e "\n"
  echo -e "Input data: " $1
  echo -e "Output: " Step1_Results/"$1"
  echo -e "Temporary workdirs: " Step1_WorkDirs/"$1"
  BASEDIR=$(pwd)
  resolve_nextflow() {
    local nf_bin
    nf_bin="$(command -v nextflow || true)"
    if [ -z "$nf_bin" ] && [ -x /usr/local/bin/nextflow ]; then
      nf_bin="/usr/local/bin/nextflow"
    fi
    if [ -z "$nf_bin" ]; then
      echo "ERROR: nextflow not found in PATH." >&2
      exit 127
    fi
    echo "$nf_bin"
  }

  run_nextflow() {
    local nf_bin
    nf_bin="$(resolve_nextflow)"
    if command -v stdbuf >/dev/null 2>&1; then
      stdbuf -oL -eL "$nf_bin" "$@"
    else
      "$nf_bin" "$@"
    fi
  }
  ## Create output directories
  mkdir -p Step1_Results/"$1"
  mkdir -p Step1_WorkDirs/"$1"
  cd Step1_WorkDirs/"$1"


## Step-1 - with pre-demultiplexed data
run_nextflow run /opt/pipelines/NextITS/main.nf \
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
  2>&1 | sed -E 's/\x1B\[[0-9;]*[A-Za-z]//g' \
       | tr -d '\r' \
       | tr -d '▒░' \
       | tee -a "$BASEDIR"/Input/Nextflow__"$1".log
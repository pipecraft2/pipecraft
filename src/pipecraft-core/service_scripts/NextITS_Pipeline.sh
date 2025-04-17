#!/bin/bash


export NXF_HOME="/input/.nextflow"
mkdir -p $NXF_HOME
BASEDIR=$(pwd)

fix_permissions() {
  # Try different possible locations for the NextITS scripts
  for dir in \
    "$NXF_HOME/assets/vmikk/NextITS/bin" \
    "/input/.nextflow/assets/vmikk/NextITS/bin" \
    "$HOME/.nextflow/assets/vmikk/NextITS/bin" \
    "./work/*/vmikk/NextITS/bin"
  do
    if [ -d "$dir" ]; then
      echo "Setting permissions for scripts in $dir"
      find "$dir" -name "*.R" -exec chmod +x {} \; 2>/dev/null
      find "$dir" -name "*.py" -exec chmod +x {} \; 2>/dev/null
      find "$dir" -name "*.sh" -exec chmod +x {} \; 2>/dev/null
      # Also fix line endings in case they're causing issues
      find "$dir" -name "*.R" -exec sed -i 's/\r$//' {} \; 2>/dev/null
      find "$dir" -name "*.py" -exec sed -i 's/\r$//' {} \; 2>/dev/null
      find "$dir" -name "*.sh" -exec sed -i 's/\r$//' {} \; 2>/dev/null
    fi
  done
}

fix_permissions

ls -la

## Run Step-1 for all sequencing runs
find /input/Input/ -type d -not -path /input/Input/ | sort \
  | parallel -j1 --joblog Step1.log \
  "/scripts/submodules/NextITS_Step1.sh {/}"




## Step-2 - standard VSEARCH clustering
nextflow run \
  vmikk/NextITS -r main \
  -resume \
  --storagemode "copy" \
  -params-file /scripts/NextFlowConfig.json \
  --step "Step2" \
  --data_path  /input/Step1_Results \
  --outdir     /input/Step2_Results \
  -work-dir    /input/Step2_WorkDir \
  --tracedir   /input/Step2_WorkDir/pipeline_info \
  -ansi-log    false \
  >> Nextflow__Step2.log

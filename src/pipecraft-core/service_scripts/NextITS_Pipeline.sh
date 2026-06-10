#!/bin/bash


export NXF_HOME="/opt/software/conda/bin"
export NXF_ANSI_LOG="false"
export NXF_LOG_COLOR="false"
export NXF_ANSI="false"
export TERM="dumb"
BASEDIR=$(pwd)



ls -la

#######################################################################
## Chimera handling
## Chimera mode is derived only from whether a chimera_db was supplied
## (there is no front-end chimera_methods parameter):
##   - no DB    -> "denovo"
##   - DB given -> "ref,denovo", and the DB is converted to VSEARCH UDB
##                 (required by NextITS) and repointed in the params file.
#######################################################################
CONFIG_FILE="/scripts/NextFlowConfig.json"

# Read a flat JSON string value from the params file
json_get() { sed -n 's/.*"'"$1"'"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$CONFIG_FILE" | head -n1; }

# Set a JSON string value: replace it if the key exists, else insert after "{"
json_set() {
  if grep -q "\"$1\"" "$CONFIG_FILE"; then
    sed -i 's|\("'"$1"'"[[:space:]]*:[[:space:]]*"\)[^"]*"|\1'"$2"'"|' "$CONFIG_FILE"
  else
    sed -i 's|^{|{"'"$1"'":"'"$2"'",|' "$CONFIG_FILE"
  fi
}

chimera_db="$(json_get chimera_db)"

if [ -z "$chimera_db" ] || [ "$chimera_db" = "undefined" ]; then
  chimera_methods="denovo"
else
  chimera_methods="ref,denovo"

  # The front-end /extraFilesN index may be wrong; locate the file by basename
  db_basename="$(basename "$chimera_db")"
  for d in /extraFiles /extraFiles[0-9]*; do
    [ -f "$d/$db_basename" ] && { chimera_db="$d/$db_basename"; break; }
  done
  [ -f "$chimera_db" ] || { echo "ERROR: chimera database not found: $db_basename" >&2; exit 1; }

  # Convert to UDB unless it already is one
  if [ "$(printf '%s' "${chimera_db##*.}" | tr '[:upper:]' '[:lower:]')" = "udb" ]; then
    final_db="$chimera_db"
  else
    # Write/reuse the UDB next to the source FASTA so it survives across runs
    final_db="${chimera_db%.*}.udb"
    if [ -s "$final_db" ]; then
      echo "# Reusing existing UDB next to FASTA: $final_db"
    else
      echo "# Converting chimera DB to UDB: $chimera_db -> $final_db"
      vsearch --makeudb_usearch "$chimera_db" --output "$final_db"
      { [ $? -eq 0 ] && [ -s "$final_db" ]; } || { echo "ERROR: failed to build UDB from $chimera_db" >&2; exit 1; }
    fi
  fi
  json_set chimera_db "$final_db"
fi

json_set chimera_methods "$chimera_methods"
echo "# chimera_methods=$chimera_methods; final config:"; cat "$CONFIG_FILE"

## Run Step-1 for all sequencing runs
find /Input/ -mindepth 1 -maxdepth 1 -type d \
  ! -name ".nextflow" \
  ! -name "Step1_Results" \
  ! -name "Step1_WorkDirs" \
  ! -name "Step2_Results" \
  ! -name "Step2_WorkDir" \
  | sort \
  | parallel -j1 --joblog /Input/Step1.log \
  "/scripts/submodules/NextITS_Step1.sh {/}"




## Step-2 - standard VSEARCH clustering

stdbuf -oL -eL \
  nextflow run /opt/pipelines/NextITS/main.nf \
  -resume \
  --storagemode "copy" \
  -params-file /scripts/NextFlowConfig.json \
  --step "Step2" \
  --data_path  /Input/Step1_Results \
  --outdir     /Input/Step2_Results \
  -work-dir    /Input/Step2_WorkDir \
  --tracedir   /Input/Step2_WorkDir/pipeline_info \
  -ansi-log    false \
  2>&1 | sed -E 's/\x1B\[[0-9;]*[A-Za-z]//g' \
       | tr -d '\r' \
       | tr -d '▒░' \
       | tee -a /Input/Nextflow__Step2.log

#!/bin/bash


export NXF_HOME="/Input/.nextflow"
mkdir -p "$NXF_HOME"
export NXF_ANSI_LOG="false"
export NXF_LOG_COLOR="false"
export NXF_ANSI="false"
export TERM="dumb"
BASEDIR=$(pwd)



ls -la

#######################################################################
## Chimera handling
## The app writes chimera_methods and the final chimera_db path into
## /scripts/NextFlowConfig.json (read-only here). The chimera_db environment
## variable is the selected file. A FASTA is converted to /dbcache/db.udb
## before Nextflow reads the params file.
#######################################################################
source_db="${chimera_db:-}"

if [ -n "$source_db" ] && [ "$source_db" != "undefined" ]; then
  [ -f "$source_db" ] || { echo "ERROR: chimera database not found: $(basename "$source_db")" >&2; exit 1; }

  if [ "$(printf '%s' "${source_db##*.}" | tr '[:upper:]' '[:lower:]')" != "udb" ]; then
    if [ ! -d /dbcache ] || [ ! -w /dbcache ]; then
      echo "ERROR: database cache is not mounted at /dbcache" >&2
      exit 1
    fi
    source /scripts/submodules/db_cache.sh
    final_db="/dbcache/db.udb"
    if db_cache_is_complete /dbcache; then
      echo "# Using cached UDB for $(basename "$source_db")"
    else
      echo "# Converting chimera DB to UDB: $(basename "$source_db") (saved for future runs)"
      db_cache_reset /dbcache
      vsearch --makeudb_usearch "$source_db" --output "$final_db"
      { [ $? -eq 0 ] && [ -s "$final_db" ]; } || { echo "ERROR: failed to build UDB from $source_db" >&2; exit 1; }
      db_cache_mark_complete /dbcache
    fi
  fi
fi

echo "# final config:"; cat /scripts/NextFlowConfig.json

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

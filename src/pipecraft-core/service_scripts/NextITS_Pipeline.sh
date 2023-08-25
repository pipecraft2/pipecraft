#!/bin/bash

## Run Step-1 for all sequencing runs
find /input/Input/ -type d -not -path /input/Input/ | sort \
  | parallel -j1 --joblog Step1.log \
  "/scripts/submodules/NextITS_Step1.sh {/}"


## Run Step 2
nextflow run /scripts/NextITS/Step2_AggregateRuns.nf \
  -resume \
  -params-file /scripts/NextFlowConfig.json \
  -profile     docker \
  --merge_replicates false \
  --data_path  /input/Step1_Results \
  --outdir     /input/Step2_Results \
  -work-dir    /input/Step2_WorkDir \
  --tracedir   /input/Step2_WorkDir/pipeline_info \
  -ansi-log    false \
  >> Nextflow__Step2.log


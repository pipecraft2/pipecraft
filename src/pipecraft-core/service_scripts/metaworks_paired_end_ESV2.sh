#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.11.2

cd /MetaWorks1.11.2
snakemake --jobs 2 --snakefile snakefile_ESV --configfile /input/config_testing_COI_data.yaml

echo "workingDir=metaworks_output"
echo "fileFormat=fastq"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"

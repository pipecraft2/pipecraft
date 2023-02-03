#!/bin/bash

echo '1'
#exec $BASH
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.11.2
echo "2"
snakemake --jobs 1 --snakefile hello_world.txt
echo "3"
echo "workingDir=metaworks_output"
echo "fileFormat=fastq"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"

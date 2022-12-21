#!/bin/sh

exec $SHELL
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.11.2
snakemake -h
snakemake --jobs 1 --snakefile /scripts/hello_world.txt
echo "workingDir=metaworks_output"
echo "fileFormat=fastq"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"
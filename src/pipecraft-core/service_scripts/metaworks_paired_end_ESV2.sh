#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.11.2

set -e 

cd /MetaWorks1.11.2

#specify output dir in the config file
#mkdir -p /input/docmw2

sed -e 's/dir:.*/dir: "\/input\/docmw2"/' config_testing_COI_data.yaml > config.yaml

#specify database
sed -i 's/\/path\/to\/mydata_trained\/rRNAClassifier.properties/\/extraFiles2\/rRNAClassifier.properties/' config.yaml 
#run snakemake
snakemake --jobs 2 --snakefile snakefile_ESV --configfile config.yaml

#echo "workingDir=metaworks_output"
echo "fileFormat=fastq"
echo "dataFormat=demultiplexed"
echo "readType=paired_end"

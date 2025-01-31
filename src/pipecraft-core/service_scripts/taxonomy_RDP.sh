#!/bin/bash

#Taxonomy annotation with RDP Classifier
#Input = single-end fasta file + database file.

################################################
###Third-party applications:
# The RDP classifier v2.13 [via metaworks container]
################################################

#load variables
regex='[^/]*$'
db1_temp=$(echo $database | grep -oP "$regex")
db1=$(printf "/extraFiles/$db1_temp")
echo "db1 = $db1"
confidence=$"--conf ${confidence}"  # default is 0.8. Assignment confidence cutoff used to determine the assignment count for each taxon. Range [0-1]
mem=$"-Xmx${mem}g"                  # default is 8GB. The amount of memory to allocate to the RDP classifier

# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/taxonomy_out"

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
#initiate conda env on container
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.12.0

### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
check_gz_zip_SE
### Check input formats (fasta supported)
check_extension_fasta
### Select last fasta file in the folder as input for RDP
for file in *.$fileFormat; do
	IN=$(echo $file)
done
echo "input = $IN"

## Perform taxonomy annotation
printf '%s\n' "Running RDP "
checkerror=$(rdp_classifier $mem classify \
			--shortseq_outfile $output_dir/shortseqs.txt \
			$confidence \
			-t $db1 \
			-o $output_dir/taxonomy.txt \
			$IN 2>&1)
check_app_error

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"

if [[ $debugger != "true" ]]; then
	if [[ -d tempdir ]];then
		rm -r tempdir
	fi
	if [[ -d tempdir2 ]];then
		rm -r tempdir2
	fi
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
db_x=$(echo $db1 | sed -e 's/\/extraFiles\///')
printf "# Taxonomy was assigned using RDP classifier (see 'Core command' below for the used settings).

Query    = $IN
Database = $db_x

# taxonomy.txt  = RDP classifier results; tab-delimited file with taxonomic ranks and associated bootstrap values
# shortseqs.txt = sequence names that are too short to be classified

Core command -> 
rdp_classifier $mem classify --shortseq_outfile shortseqs.txt $confidence -t $db_x -o taxonomy.txt $IN

Total run time was $runtime sec.

################################################
###Third-party applications:
    #The RDP classifier v2.13:
	#citation: Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology. 73(16), 5261-5267. doi:10.1128/AEM.00062-07
#####################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

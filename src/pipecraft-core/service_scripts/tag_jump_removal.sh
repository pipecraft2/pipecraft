#!/bin/bash

#Using UNCROSS2 to filter out putative tag-jumps.
# Input = OTU_table

##################################
###Third-party applications:
#R data.table, ggplot2
##################################

# $OTU_table = input OTU_table
regex='[^/]*$'
OTU_table_file_path=$(echo $OTU_table | grep -oP "$regex")
OTU_table_temp=$(basename $OTU_table_file_path) #basename, needed for macOS
OTU_table_file=$(printf "/extraFiles/$OTU_table_temp")

# f_value = 
# p_value = 

#output dir
output_dir=$"/input/"
# Source for functions
source /scripts/submodules/framework.functions.sh

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
printf "# Running tag-jumps filtering \n"
Rlog=$(Rscript /scripts/submodules/tag_jump_removal.R $OTU_table_file $f_value $p_value 2>&1)
echo $Rlog > $output_dir/tag-jumps_filt.log 
wait
printf "\n tag-jumps filtering completed \n"
#format R-log file
sed -i "s/;; /\n/g" $output_dir/tag-jumps_filt.log 

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
if [[ $debugger != "true" ]]; then
    rm $output_dir/TagJump_OTUs.RData 
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
in=$(echo $OTU_table_file | sed -e 's/\/extraFiles\///')

printf "# Used UNCROSS2 algorithm with f-parameter = $f_value and p-parameter = $p_value
    f-parameter defines the expected cross-talk rate. Default is 0.01 (equivalent to 1%). A higher value enforces stricter filtering.
    p-parameter controls the severity of tag-jump removal. It adjusts the exponent in the UNCROSS formula. Default is 1. Opt for e.g. 0.5 or 0.3 to steepen the curve.

Input feature (OTU/ASV) table = $in.

Generated files:
# TagJumpFiltered_FeatureTable = output table where tag-jupms have been filtered out
# TagJump_plot.pdf             = 
# TagJump_stats.txt            = statistics about 
                                Total_reads             : number of reads in the input table 
                                Number_of_TagJump_Events: number of cases where tag-jumps were detected
                                TagJump_reads           : number of reads detected as tag-jumps
                                ReadPercent_removed     : percentence of the reads removed 

Total run time was $runtime sec.

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#UNCROSS2
    citation: R.C. Edgar (2018), UNCROSS2: identification of cross-talk in 16S rRNA OTU tables, https://doi.org/10.1101/400762
##################################################################" > $output_dir/README_TagJumps.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
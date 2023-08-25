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

# f_value       = 
# p_value       = 

#output dir
output_dir=$"/input/"
# Source for functions
source /scripts/submodules/framework.functions.sh

#############################
### Start of the workflow ###
#############################
printf "# Running tag-jumps filtering \n"
Rlog=$(Rscript /scripts/submodules/tag_jump_removal.R $OTU_table_file $f_value $p_value 2>&1)
echo $Rlog > $output_dir/tag-jumps_filt.log 
wait
printf "\n tag-jumps filtering completed \n"
#format R-log file
sed -i "s/;; /\n/g" $output_dir/tag-jumps_filt.log 


#rm $output_dir/TagJump_OTUs.RData 

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
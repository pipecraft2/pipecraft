#!/bin/bash

# Collapse no mismatch ASVs and/or filter ASVs by length in a dada2 made ASV table 
    #Script for the DADA2 ASVs workflow
# Input = DADA2 ASV_table file in rds format (saved by the DADA2 ASVs workflow)

##########################################################
###Third-party applications:
#dada2 v1.20
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
##########################################################
set -e
#load variables
collapseNoMismatch=${collapseNoMismatch}
minOverlap=${minOverlap}
vec=${vec}
by_length=${by_length}

#Source for functions
source /scripts/submodules/framework.functions.sh

#Quit if collapseNoMismatch = false, and by_length filtering value is 0.
if [[ $collapseNoMismatch != "true" ]]; then
    if [[ $by_length == 0 ]]; then
        printf "collapseNoMismatch = false, and by_length filtering value is 0. Nothing to be done. Quitting\n"
        end_process
    fi
fi

#output dir
output_dir=$"/input/filtered_table"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir -p $output_dir

#############################
### Start of the workflow ###
#############################
start=$(date +%s)

### Filtering the ASV table
printf "Filtering the ASV table ... \n"
Rlog=$(Rscript /scripts/submodules/dada2_table_filtering.R 2>&1)
echo $Rlog > $output_dir/R_run.log 
wait

# Count ASVs
ASVs_collapsed=$(grep -c "^>" $output_dir/ASVs_collapsed.fasta)
if [[ $ASVs_collapsed == "" ]]; then
    ASVs_collapsed=$"0"
fi

ASVs_lenFilt=$(grep -c "^>" $output_dir/ASVs_lenFilt.fasta)
if [[ -f $output_dir/a.txt ]]; then
    ASVs_lenFilt_result_table=$"None of the ASVs were filtered out based on the length filter ($by_length bp); no files generated"
    rm $output_dir/a.txt
elif [[ $ASVs_lenFilt == "" ]]; then
    ASVs_lenFilt_result_table=$"All ASVs were filtered out based on the length filter ($by_length bp); no files generated"
elif [[ $ASVs_lenFilt != "" ]]; then
    ASVs_lenFilt_result_table=$"ASV table after discarding < $by_length bp ASVs. Contains $ASVs_lenFilt ASVs"
    ASVs_lenFilt_result_fasta=$"ASV sequences after discarding < $by_length bp ASVs. Contains $ASVs_lenFilt ASVs"
fi

# Add "ASV" as a 1st col name
if [[ -s $output_dir/ASV_table_lenFilt.txt ]]; then
    sed -i "1 s|^|ASV|" $output_dir/ASV_table_lenFilt.txt
fi
if [[ -s $output_dir/ASVs_table_collapsed.txt ]]; then
    sed -i "1 s|^|ASV|" $output_dir/ASVs_table_collapsed.txt
fi

#rm 
if [[ -f $output_dir/R_run.log ]]; then
    rm -f $output_dir/R_run.log
fi

#Make README.txt file
end=$(date +%s)
runtime=$((end-start))
if [[ $collapseNoMismatch == "true" ]]; then
    printf "# Filtering the of the dada2 ASV table.
    Note that ASV headers between tables/fastas do not match. Named from ASV_1 to ASV_N according to ASV lis

    Files in 'filtered_table':
    # ASVs_table_collapsed.txt = ASV table after collapsing identical ASVs. Contains $ASVs_collapsed ASVs
    # ASVs_collapsed.fasta     = ASV sequences after collapsing identical ASVs. Contains $ASVs_collapsed ASVs
    # ASV_table_collapsed.rds  = ASV table in RDS format after collapsing identical ASVs. 

    # If length filtering was applied (if 'by length' setting > 0) [performed after collapsing identical ASVs]:
    # ASV_table_lenFilt.txt    = $ASVs_lenFilt_result_table
    # ASVs_lenFilt.fasta       = $ASVs_lenFilt_result_fasta


    Total run time was $runtime sec.
    ##################################################################
    ###Third-party applications for this process [PLEASE CITE]:
    #dada2 v1.20
        #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
        #https://github.com/benjjneb/dada2
    ########################################################" > $output_dir/README.txt
else
    printf "# Filtering the of the dada2 ASV table.
    Note that ASV headers between tables/fastas do not match. Named from ASV_1 to ASV_N according to ASV lis

    Files in 'filtered_table':
    # ASV_table_lenFilt.txt    = $ASVs_lenFilt_result_table
    # ASVs_lenFilt.fasta       = $ASVs_lenFilt_result_fasta

    Total run time was $runtime sec.
    ##################################################################
    ###Third-party applications for this process [PLEASE CITE]:
    #dada2 v1.20
        #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
        #https://github.com/benjjneb/dada2
    ########################################################" > $output_dir/README.txt
fi

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Check README_ASVtabFilt.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "dataFormat=$dataFormat"
echo "readType=single_end"
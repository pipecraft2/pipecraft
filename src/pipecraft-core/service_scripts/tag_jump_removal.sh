#!/bin/bash

#Using UNCROSS2 to filter out putative tag-jumps.
# Input = OTU_table; TAB-DELIMITED OTU/ASV table, 
                # where the 1st column is the OTU/ASV IDs and the following columns represent samples; 
                # 2nd column may be Sequence column, with the colName 'Sequence'
# edit 21.01.2025: integration to full pipelines and multiRunDir handling

##############################
### Third-party applications:
# R data.table, ggplot2
##############################

# Source for functions
source /scripts/submodules/framework.functions.sh

# input from Quick Tools (standalone)
    # below is the input from full pipelines
if [[ -n $OTU_table ]] && [[ -n $f_value ]] && [[ -n $p_value ]]; then
    regex='[^/]*$'
    OTU_table_file_path=$(echo $OTU_table | grep -oP "$regex")
    OTU_table_temp=$(basename $OTU_table_file_path)
    OTU_table_file=$(printf "/extraFiles/$OTU_table_temp")
    printf "\n running UNCROSS2 as a standalone process, input = $OTU_table_file\n"
    #output dir
    output_dir="/input/"
    # Validate f_value is a valid number between 0 and 1
    if ! [[ $f_value =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$f_value < 0" | bc -l) )) || (( $(echo "$f_value > 10" | bc -l) )); then
        log_error "f_value must be a number between 0 and 10. Received value: $f_value
        >Quitting"
        end_process
    fi
    # Validate p_value is a positive number
    if ! [[ $p_value =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$p_value <= 0" | bc -l) )); then
        log_error "p_value must be a positive number. Received value: $p_value
        >Quitting"
        end_process
    fi
fi 

# R dependency check
check_r_dependencies() {
    Rscript -e 'if(!require("data.table")) quit(status=1); if(!require("ggplot2")) quit(status=1)' 2>/dev/null
    if [ $? -ne 0 ]; then
        log_error "Required R packages (data.table and/or ggplot2) are not installed.
        Please install missing R dependencies.
        >Quitting"
        end_process
    fi
}

# Check if I need to work with multiple or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "Tag-jumps filtering with UNCROSS2 for multiple sequencing runs in multiRunDir"
    echo "Process = tag-jumps filtering"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | sed -e "s/^\.\///")
    echo "Working in dirs:"
    echo $DIRS
    multiDir="TRUE"
    export multiDir
else
    cd /input
    echo "Working with individual sequencing run"
    echo "Process = tag-jumps filtering"
    DIRS="/input"
    printf "\n workingDirs = $DIRS \n"
    multiDir="FALSE"
    export multiDir
fi

#############################
### Start of the workflow ###
#############################
### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start=$(date +%s)
    cd $seqrun

    # Multi-sequencing run (full pipeline)
    if [[ $multiDir == "TRUE" ]]; then
        # Check for input table; define output_dir and workingDir
        if [[ -f "/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/ASVs_table.txt" ]]; then
            OTU_table_file="/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2/ASVs_table.txt"
            output_dir="/input/multiRunDir/${seqrun%%/*}/ASVs_out.dada2"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with ASV table, input = $OTU_table_file\n"
        elif [[ -f "/input/multiRunDir/${seqrun%%/*}/clustering_out/OTU_table.txt" ]]; then
            OTU_table_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/OTU_table.txt"
            output_dir="/input/multiRunDir/${seqrun%%/*}/clustering_out/"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with OTU table, input = $OTU_table_file\n"
        elif [[ -f "/input/multiRunDir/${seqrun%%/*}/clustering_out/zOTU_table.txt" ]]; then
            OTU_table_file="/input/multiRunDir/${seqrun%%/*}/clustering_out/zOTU_table.txt"
            output_dir="/input/multiRunDir/${seqrun%%/*}/clustering_out/"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with zOTU table, input = $OTU_table_file\n"
        else
            printf '%s\n' "ERROR]: Could not find input table.
            Looked for:
            - /ASVs_out.dada2/ASVs_table.txt
            - /clustering_out/OTU_table.txt
            - /clustering_out/zOTU_table.txt
            Please check if ASV/OTU table exists.
            >Quitting" >&2
            end_process
        fi
    # Single sequencing run (full pipeline)
    elif [[ $multiDir == "FALSE" ]]; then   
        # Check for input table; define output_dir and workingDir
        if [[ -f "/input/ASVs_out.dada2/ASVs_table.txt" ]]; then
            OTU_table_file="/input/ASVs_out.dada2/ASVs_table.txt"
            output_dir="/input/ASVs_out.dada2"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with ASV table, input = $OTU_table_file\n"
        elif [[ -f "/input/clustering_out/OTU_table.txt" ]]; then
            OTU_table_file="/input/clustering_out/OTU_table.txt"
            output_dir="/input/clustering_out/"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with OTU table, input = $OTU_table_file\n"
        elif [[ -f "/input/clustering_out/zOTU_table.txt" ]]; then
            OTU_table_file="/input/clustering_out/zOTU_table.txt"
            output_dir="/input/clustering_out/"
            export output_dir
            workingDir=$output_dir
            export workingDir
            printf "\n running UNCROSS2 with zOTU table, input = $OTU_table_file\n"
        else
            printf '%s\n' "ERROR]: Could not find input table.
            Looked for:
            - /ASVs_out.dada2/ASVs_table.txt
            - /clustering_out/OTU_table.txt
            - /clustering_out/zOTU_table.txt
            Please check if ASV/OTU table exists.
            >Quitting" >&2
            end_process
        fi
    fi
 
    ### Process samples with UNCROSS2 in R
    printf "# Running UNCROSS2 for $seqrun\n "
    Rlog=$(Rscript /scripts/submodules/tag_jump_removal.R $OTU_table_file $f_value $p_value 2>&1)
    # Check if R script executed successfully
    if [ $? -ne 0 ]; then
        log_error "tag-jumps filtering R script failed with the following error:
        $Rlog
        Please check the parameters and input file.
        >Quitting"
        end_process
    fi
    echo $Rlog > $output_dir/tag-jumps_filt.log 
    wait
    
    # Check if output files were created
    if [ ! -f "$output_dir/TagJumpFiltered_FeatureTable.txt" ]; then
        log_error "UNCROSS2 did not generate the expected output file.
        Please check the log file at $output_dir/tag-jumps_filt.log
        >Quitting"
        end_process
    fi

    printf "\n tag-jumps filtering completed \n"
    # format R-log file
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
        f-parameter defines the expected cross-talk rate. Default is 0.01 (equivalent to 1%%). A higher value enforces stricter filtering.
        p-parameter controls the severity of tag-jump removal. It adjusts the exponent in the UNCROSS formula. Default is 1. Opt for e.g. 0.5 or 0.3 to steepen the curve.

    Input feature (OTU/ASV) table = $in.

    Generated files:
    # TagJumpFiltered_FeatureTable = output table where tag-jupms have been filtered out
    # TagJump_plot.pdf             = illustration about the presence of tag-jumps based on the selected parameters
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

    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    else
        cd /input
    fi
done

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
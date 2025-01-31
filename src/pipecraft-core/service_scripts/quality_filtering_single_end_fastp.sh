#!/bin/bash

# Quality filter SINGLE-END sequencing data with fastp
# Input = single-end fastq files

################################################
###Third-party applications:
#fastp v0.23.2
    #citation: Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
    #Copyright (c) 2016 OpenGene - Open Source Genetics Toolbox
    #Distributed under the MIT License
    #https://github.com/OpenGene/fastp
#pigz v2.4
################################################

#load variables
window_size=$"--cut_window_size ${window_size}"
required_qual=$"--cut_mean_quality ${required_qual}"
min_qual=$"--qualified_quality_phred ${min_qual}"
min_qual_thresh=$"--unqualified_percent_limit ${min_qual_thresh}"
maxNs=$"--n_base_limit ${maxNs}"
min_length=$"--length_required ${min_length}"
max_length=$"--length_limit ${max_length}"
trunc_length=$"--max_len1 ${trunc_length}"
aver_qual=$"--average_qual ${aver_qual}"
cores=$"--thread ${cores}"

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/qualFiltered_out"

#additional options, if selection != undefined
low_complex_filt=$low_complexity_filter
if [[ $low_complex_filt == null ]] || [[ -z $low_complex_filt ]]; then
    low_complexity_filter=$""
else
    low_complexity_filter=$"--low_complexity_filter --complexity_threshold $low_complex_filt"
fi

trim_polyG=${trim_polyG}
if [[ $trim_polyG == null ]] || [[ -z $trim_polyG ]]; then
    trim_polyG=$"--disable_trim_poly_g "
else
    trim_polyG=$"--trim_poly_g --poly_g_min_len $trim_polyG"
fi

trim_polyX=${trim_polyX}
if [[ $trim_polyX == null ]] || [[ -z $trim_polyX ]]; then
    trim_polyX=$""
else
    trim_polyX=$"--trim_poly_x --poly_x_min_len $trim_polyX"
fi

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env
### Process samples
for file in *.$fileFormat; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
    ## Preparing files for the process
    printf "\n____________________________________\n"
    printf "Processing $input ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_gz_zip_SE
    ### Check input formats (fastq supported)
    check_extension_fastq

    ###############################
    ### Start quality filtering ###
    ###############################
    checkerror=$(fastp --in1 $input.$extension \
                       --out1 $output_dir/$input.$extension \
                       $window_size \
                       $required_qual \
                       $min_qual \
                       $min_qual_thresh \
                       $trim_polyG \
                       $trim_polyX \
                       $maxNs \
                       $min_length \
                       $max_length \
                       $trunc_length \
                       $aver_qual \
                       $cores \
                       --html $output_dir/fastp_report/$input.html \
                       --disable_adapter_trimming \
                       $low_complexity_filter 2>&1)
                       check_app_error
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering was performed using fastp (see 'Core command' below for the used settings).

Files in 'qualFiltered_out':
----------------------------
# *.fastq               = quality filtered sequences per sample.
# seq_count_summary.txt = summary of sequence counts per sample.

Core command -> 
fastp --in1 input --out1 output $window_size $required_qual $min_qual $min_qual_thresh $trim_polyG $trim_polyX $maxNs $min_length $max_length $trunc_length $aver_qual $cores --html fastp_report/sample_name.html --disable_adapter_trimming $low_complexity_filter

Total run time was $runtime sec.

##############################################
###Third-party applications for this process:
#fastp v0.23.2
    #citation: Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
    #https://github.com/OpenGene/fastp
##############################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

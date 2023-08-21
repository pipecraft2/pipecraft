#!/bin/bash

# Quality filter PAIRED-END sequencing data with fastp
# Input = paired-end fastq files

##########################################################
###Third-party applications:
#fastp v0.23.2
    #citation: Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
    #Copyright (c) 2016 OpenGene - Open Source Genetics Toolbox
    #Distributed under the MIT License
    #https://github.com/OpenGene/fastp
#pigz v2.4
##########################################################

#load variables
extension=$fileFormat
window_size=$"--cut_window_size ${window_size}"
required_qual=$"--cut_mean_quality ${required_qual}"
min_qual=$"--qualified_quality_phred ${min_qual}"
min_qual_thresh=$"--unqualified_percent_limit ${min_qual_thresh}"
maxNs=$"--n_base_limit ${maxNs}"
min_length=$"--length_required ${min_length}"

max_length=$"--length_limit ${max_length}"
trunc_length_R1=$"--max_len1 ${trunc_length}"
trunc_length_R2=$"--max_len2 ${trunc_length}"
aver_qual=$"--average_qual ${aver_qual}"
cores=$"--thread ${cores}"

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/qualFiltered_out"

#additional options, if selection != undefined
low_complex_filt=${low_complexity_filter}
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
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env
#Make dir for fastp_reports
mkdir -p $output_dir/fastp_report
### Process samples
while read LINE; do
    #Read in R1 and R2 file names; without extension
    inputR1=$(echo $LINE | sed -e "s/.$extension//")
    inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
    sample_name=$(echo $inputR1 | sed -e 's/R1/fastp_report/')
    ## Preparing files for the process
    printf "\n____________________________________\n"
    printf "Processing $inputR1 and $inputR2 ...\n"
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        
    check_gz_zip_PE
    ### Check input formats (fastq supported)
    check_extension_fastq

    ###############################
    ### Start quality filtering ###
    ###############################
    checkerror=$(fastp --in1 $inputR1.$extension \
                       --in2 $inputR2.$extension \
                       --out1 $output_dir/$inputR1.$extension \
                       --out2 $output_dir/$inputR2.$extension \
                       $window_size \
                       $required_qual \
                       $min_qual \
                       $min_qual_thresh \
                       $maxNs \
                       $min_length \
                       $max_length \
                       $trunc_length_R1 \
                       $trunc_length_R2 \
                       $aver_qual \
                       $trim_polyG \
                       $trim_polyX \
                       $cores \
                       --html $output_dir/fastp_report/$sample_name.html \
                       --disable_adapter_trimming \
                       $low_complexity_filter 2>&1)
                       check_app_error
done < tempdir2/paired_end_files.txt

#remove json reports
if [[ -s fastp.json ]]; then
    rm fastp.json
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
clean_and_make_stats
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# Quality filtering with fastp.

Files in 'qualFiltered_out':
# *.fastq               = quality filtered sequences per sample.
# seq_count_summary.txt = summary of sequence counts per sample.

Core command -> 
fastp --in1 inputR1 --in2 inputR2 --out1 outputR1 --out2 outputR2 $window_size $required_qual $min_qual $min_qual_thresh $trim_polyG $trim_polyX $maxNs $min_length $max_length $trunc_length_R1 $trunc_length_R2 $aver_qual $cores --html fastp_report/sample_name.html --disable_adapter_trimming $low_complexity_filter

Total run time was $runtime sec.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#fastp v0.23.2
    #citation: Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
    #https://github.com/OpenGene/fastp
########################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"

echo "readType=paired_end"

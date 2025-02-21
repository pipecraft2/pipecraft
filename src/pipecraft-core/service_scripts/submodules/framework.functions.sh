#!/bin/bash

# Set of functions for PipeCraft (2.0) workflows, for checking data integrity.
##############################
### Error logging function ###
##############################
log_error() {
    local error_message="$1"
    local timestamp=$(date '+%Y-%m-%d_%H.%M.%S')
    local log_file="$output_dir/log_$timestamp.log"
      
    # Log to file
    echo "[$timestamp] ERROR: $error_message" >> "$log_file"
    # Print to stderr
    printf '%s\n' "ERROR]: $error_message" >&2
}

###############################
### Quit process upon ERROR ###
###############################
function end_process () {
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir" ]]; then
        rm -rf tempdir
    fi
    if [[ -d "tempdir2" ]]; then
        rm -rf tempdir2
    fi
fi
exit 1
}
export -f end_process

#######################################
### Check if APP run was successful ###
#######################################
function check_app_error () {
if [[ "$?" = "0" ]]; then
    :
else
    printf '%s\n' "ERROR]: $checkerror" >&2
    end_process
fi
}
export -f check_app_error

##################################################
### Function to convert IUPAC codes in primers ###
##################################################
function convert_IUPAC () {
echo $1 | \
if grep -q -E "R|Y|S|W|K|M|B|D|H|V|N|I"; then
    #define IUPAC codes
    R=$"[AG]"
    Y=$"[CT]"
    S=$"[GC]"
    W=$"[AT]"
    K=$"[GT]"
    M=$"[AC]"
    B=$"[CGT]"
    D=$"[AGT]"
    H=$"[ACT]"
    V=$"[ACG]"
    N=$"[ATGC]"
    I=$"[ATGC]"
    #replace IUPAC codes
    primer=$(echo $1 | \
    sed -e "s/R/$R/g; s/Y/$Y/g; \
    s/S/$S/g; s/W/$W/g; s/K/$K/g; \
    s/M/$M/g; s/B/$B/g; s/D/$D/g; \
    s/H/$H/g; s/V/$V/g; s/N/$N/g; \
    s/I/$I/g")
    #return convered primer
    echo $primer
else
    #return original primer when no IUPAC codes were detected
    echo $1
fi
}

#######################################################
### Check if files with specified extension exist in the dir ###
#######################################################
function first_file_check () {
count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
if (( $count != 0 )); then 
    :
else
    printf '%s\n' "ERROR]: cannot find files with specified extension '$fileFormat'
Please check the extension of your files and specify again.
>Quitting" >&2
    end_process
fi 
}

function first_file_check_clustering () {
count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
if (( $count != 0 )); then 
    :
else
    printf '%s\n' "ERROR]: cannot find input ($fileFormat files).
Please check and specify again.
>Quitting" >&2
    end_process
fi 
}


####################################################
### Check PAIRED-END data and pepare working env ###
####################################################
function prepare_PE_env () {
echo "output_dir = $output_dir"
#Remove 'old' output_dir if exist and make new empty one
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir -p $output_dir
#Make tempdir2, for seq count statistics
if [[ -d "tempdir2" ]]; then
    rm -rf tempdir2
fi
mkdir -p tempdir2
#Make a file where to read R1 and R2 file names for paired-end read processing.
touch tempdir2/files_in_folder.txt
for file in *.$fileFormat; do
    echo $file >> tempdir2/files_in_folder.txt
done
#Check for empty spaces in the files names. Replace with _
while read file; do
    if [[ $file == *" "* ]]; then
        printf '%s\n' "WARNING]: File $file name contains empty spaces, replaced 'space' with '_'" >&2
        mv "$file" "${file// /_}"
    fi
done < tempdir2/files_in_folder.txt
#Fix also names in files_in_folder file if they contained space
sed -i 's/ /_/g' tempdir2/files_in_folder.txt

#Grep R1 string is in the file names
grep "R1" < tempdir2/files_in_folder.txt > tempdir2/paired_end_files.txt || true

#Detect R1 pattern from the first R1 file
read_R1=$(head -n1 tempdir2/paired_end_files.txt | grep -o '[._]R1' | tail -n1)
if [[ -z "$read_R1" ]]; then
    printf '%s\n' "ERROR]: no paired-end read files found.
File names must contain '_R1' or '.R1' strings! (e.g. s01_R1.fastq and s01_R2.fastq)
>Quitting" >&2
    end_process
else 
    # Verify all files have the same R1 pattern
    while read file; do
        # Check for multiple R1/R2 occurrences first
        r1_count=$(echo "$file" | grep -o "R1" | wc -l)
        if [[ $r1_count -gt 1 ]]; then
            printf '%s\n' "ERROR]: File '$file' contains multiple R1 strings (i.e. there was other R1 stings besides just read identifier).
Sample names cannot contain R1 strings (e.g. 'sampleR1_R1.fastq' or 'sampleR11.R1.fastq' are not allowed).
>Quitting" >&2
            end_process
        fi

        # Then check for consistent R1 pattern
        if ! echo "$file" | grep -q "$(echo $read_R1 | sed 's/\./\\./g')"; then
            printf '%s\n' "ERROR]: File '$file' has inconsistent R1 pattern (did not contain '${read_R1}').
All files must use the same pattern ('_R1' or '.R1'). 
>Quitting" >&2
            end_process
        fi
    done < tempdir2/paired_end_files.txt
    
    echo "Detected R1 pattern: ${read_R1}"
    read_R1="${read_R1}" 
    export read_R1
fi
}

####################################################
### Check SINGLE-END data and pepare working env ###
####################################################
function prepare_SE_env () {
echo "output_dir = $output_dir"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir
if [[ -d "tempdir2" ]]; then
    rm -rf tempdir2
fi
mkdir -p tempdir2
#Make a file where to read file names for single-end read processing.
touch tempdir2/files_in_folder.txt
for file in *.$fileFormat; do
    echo $file >> tempdir2/files_in_folder.txt
done
#Check for empty spaces in the files names. Replace with _
while read file; do
    if [[ $file == *" "* ]]; then
        printf '%s\n' "WARNING]: File $file name contains empty spaces, replaced 'space' with '_'" >&2
        mv "$file" "${file// /_}"
    fi
done < tempdir2/files_in_folder.txt
# replace empty spaces in file names with _
sed -i 's/ /_/g' tempdir2/files_in_folder.txt
}

############################################
### Prepare primers file for CUT PRIMERS ###
############################################
function prepare_primers () {
# Primer arrays
fwd_primer_array=$(echo $fwd_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
rev_primer_array=$(echo $rev_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
# Forward primer(s) to fasta file
i=1
for primer in $fwd_primer_array; do
    echo ">fwd_primer$i" >> tempdir2/fwd_primer.fasta
    echo $primer >> tempdir2/fwd_primer.fasta
    ((i=i+1))
done
# Reverse primer(s) to fasta file
i=1
for primer in $rev_primer_array; do
    echo ">rev_primer$i" >> tempdir2/rev_primer.fasta
    echo $primer >> tempdir2/rev_primer.fasta
    ((i=i+1))
done
# Reverse complement FWD and REV primers
checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/fwd_primer.fasta >> tempdir2/fwd_primer_RC.fasta 2>&1)
check_app_error
checkerror=$(seqkit seq --quiet -t dna -r -p tempdir2/rev_primer.fasta >> tempdir2/rev_primer_RC.fasta 2>&1)
check_app_error

# Make linked primer files
i=1
while read LINE; do
    fwd_primer=$(echo $LINE | grep -v ">")
    if [ -z "$fwd_primer" ]; then
        :
    else
        while read LINE; do
            rev_primer_RC=$(echo $LINE | grep -v ">")
            if [ -z "$rev_primer_RC" ]; then
                :
            else
                echo ">primer$i" >> tempdir2/liked_fwd_revRC.fasta
                echo "$fwd_primer;$required_optional...$rev_primer_RC;$required_optional" >> tempdir2/liked_fwd_revRC.fasta
                ((i=i+1))
            fi
        done < tempdir2/rev_primer_RC.fasta
    fi
done < tempdir2/fwd_primer.fasta
i=1
while read LINE; do
    rev_primer=$(echo $LINE | grep -v ">")
    if [ -z "$rev_primer" ]; then
        :
    else
        while read LINE; do
            fwd_primer_RC=$(echo $LINE | grep -v ">")
            if [ -z "$fwd_primer_RC" ]; then
                :
            else
                echo ">primer$i" >> tempdir2/liked_rev_fwdRC.fasta
                echo "$rev_primer;$required_optional...$fwd_primer_RC;$required_optional" >> tempdir2/liked_rev_fwdRC.fasta
                ((i=i+1))
            fi
        done < tempdir2/fwd_primer_RC.fasta
    fi
done < tempdir2/rev_primer.fasta
}

###################################################
### Check if single-end files are compressed (decompress and check) ###
###################################################
#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
function check_gz_zip_SE () {
    #Read user specified input extension
    #If compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_compress=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $NF}';))
    if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
        pigz --decompress --force --keep $input.$fileFormat
        #Check errors
        if [[ "$?" != "0" ]]; then
            printf '%s\n' "ERROR]: $input.$fileFormat decompressing failed! File not compressed as gz or zip.
Decompressing other formats is not supported, please decompress manually.
>Quitting" >&2 
            end_process      
        fi
        extension=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $(NF-1)}';))
        export extension
        export check_compress
    else 
        extension=$fileFormat
        export extension
        export check_compress
    fi
}


###################################################
### Check if paired-end files are compressed (decompress and check) ###
###################################################
#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
function check_gz_zip_PE () {
check_compress=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $NF}';))
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
    pigz --decompress --force --keep $inputR1.$fileFormat
    #Check errors
    if [[ "$?" != "0" ]]; then
        printf '%s\n' "ERROR]: $inputR1.$fileFormat decompressing failed! File not compressed as gz or zip.
Decompressing other formats is not supported, please decompress manually.
>Quitting" >&2
        end_process
    fi
    pigz --decompress --force --keep $inputR2.$fileFormat
    #Check errors
    if [[ "$?" != "0" ]]; then
        printf '%s\n' "ERROR]: $inputR2.$fileFormat decompressing failed! File not compressed as gz or zip.
Decompressing other formats is not supported, please decompress manually.
>Quitting" >&2
        end_process
    fi
    extension=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $(NF-1)}';))
    export extension
    export check_compress
else 
    extension=$fileFormat
    export extension
    export check_compress
fi
}


#######################################
### Check file formatting for FASTQ ###
#######################################
function check_extension_fastq () {
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    :
else
    printf '%s\n' "ERROR]: $file formatting not supported here!
Supported extensions: fastq, fq (and gz or zip compressed formats).
>Quitting" >&2
    end_process
fi
}

#######################################
### Check file formatting for FASTA ###
#######################################
function check_extension_fasta () {
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    :
else
    printf '%s\n' "ERROR]: $file formatting not supported here!
Supported extensions: fasta, fas, fa (and gz or zip compressed formats).
>Quitting" >&2
    end_process
fi
}

#########################################
### Check file formatting for FASTQ/A ###
#########################################
function check_extension_fastx () {
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]] || [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported here!
Supported extensions: fastq, fq, fasta, fas, fa (and gz or zip compressed formats).
>Quitting" >&2
    end_process
fi
}

######################################################
### Check file formatting for FASTQ/A and their GZ ###
######################################################
function check_extension_fastxgz () {
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]] || [[ $extension == "fastq" ]] || [[ $extension == "fq" ]] || [[ $extension == "fasta.gz" ]] || [[ $extension == "fa.gz" ]] || [[ $extension == "fas.gz" ]] || [[ $extension == "fastq.gz" ]] || [[ $extension == "fq.gz" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported here!
Supported extensions: fastq, fq, fasta, fas, fa (and gz compressed formats).
>Quitting" >&2
    end_process
fi
}

###################################################
### Cleaning up and compiling final stats file, only for demux data ###
###################################################
function clean_and_make_stats_demux () {
#Delete empty output files
find $output_dir -empty -type f -delete
# Count input reads
touch tempdir2/seq_count.txt
if [[ -f tempdir2/paired_end_files.txt ]]; then
    while read LINE; do
        seqkit stats --threads 6 -T $LINE | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $1,$4}' >> tempdir2/seq_count.txt
    done < tempdir2/paired_end_files.txt
else
    while read LINE; do
        seqkit stats --threads 6 -T $LINE | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $1,$4}' >> tempdir2/seq_count.txt
    done < tempdir2/files_in_folder.txt
fi

### Count reads after the process
touch tempdir2/seq_count_after.txt
outfile_check=$(ls $output_dir/*.$fileFormat 2>/dev/null | wc -l)
if (( $outfile_check != 0 )); then
    if [[ -f tempdir2/paired_end_files.txt ]]; then
        for file in $output_dir/*R1.$fileFormat; do
            seqkit stats --threads 6 -T $file | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $1,$4}' | sed -e 's/demultiplex_out\///' >> tempdir2/seq_count_after.txt
            #sum all demux seqs (column 2)
            demux_sum=$(awk -F'\t' '{sum+=$2;} END{print sum;}' tempdir2/seq_count_after.txt)
        done
    else
        for file in $output_dir/*.$fileFormat; do
            seqkit stats --threads 6 -T $file | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $1,$4}' | sed -e 's/demultiplex_out\///' >> tempdir2/seq_count_after.txt
            #sum all demux seqs (column 2)
            demux_sum=$(awk -F'\t' '{sum+=$2;} END{print sum;}' tempdir2/seq_count_after.txt)
        done
    fi
else
    printf '%s\n' "ERROR]: no output files generated ($output_dir). Check index file and settings." >&2
    end_process
fi


### Compile a track reads summary file
printf "Input file:\n" > $output_dir/seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    printf "$file1\t$count1\n" >> $output_dir/seq_count_summary.txt    
done < tempdir2/seq_count.txt
printf "\nSUM of demultiplexed sequences (including 'unknown.$fileFormat'; see below)\t$demux_sum" >> $output_dir/seq_count_summary.txt 
printf "\n\nSamples\tNumber_of_seqs\n" >> $output_dir/seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    printf "$file1\t$count1\n" >> $output_dir/seq_count_summary.txt    
done < tempdir2/seq_count_after.txt

#remove "/input/" from $output_dir/seq_count_summary.txt
sed -i 's/\/input\///' $output_dir/seq_count_summary.txt

if [[ -f tempdir2/paired_end_files.txt ]]; then
    printf "\n[paired R2 file has the same number of sequencs as corresponding R1 file]\n" >> $output_dir/seq_count_summary.txt
fi 

#remove tempdir2
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir2" ]]; then
        rm -rf tempdir2
    fi
fi
}

###############################################
### Cleaning up and compiling final stats file, only for assemble PE data ###
###############################################
function clean_and_make_stats_assemble () {
#Delete empty output files
find $output_dir -empty -type f -delete
# Count input reads
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]] || [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    touch tempdir2/seq_count.txt
    while read LINE; do
        seqkit stats --threads 6 -T $LINE | awk -F'\t' 'BEGIN{OFS="\t";} FNR == 2 {print $1,$4}' >> tempdir2/seq_count.txt
    done < tempdir2/paired_end_files.txt
fi

### Count reads after the process
output_dir_for_sed=$(echo $output_dir | sed -e "s/\//\\\\\//g")
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]] || [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    touch tempdir2/seq_count_after.txt
    outfile_check=$(ls $output_dir/*.$extension 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then 
        seqkit stats --threads 6 -T $output_dir/*.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' | sed -e "s/$output_dir_for_sed//" >> tempdir2/seq_count_after.txt
    else 
        printf '%s\n' "ERROR]: no output files generated ($output_dir). Adjust settings." >&2
        end_process
    fi
fi

### Compile a track reads summary file (seq_count_summary.txt)
printf "File\tReads\tAssembled_reads\n" > $output_dir/seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}' | sed -e "s/$read_R1.*\.$fileFormat/\.$extension/")
    count1=$(echo $LINE | awk '{print $2}')
    file2=$(grep "$file1" tempdir2/seq_count_after.txt | awk '{print $1}' | awk 'BEGIN{FS="/"}{print $NF}')
    count2=$(grep "$file1" tempdir2/seq_count_after.txt | awk '{print $2}')
    if [[ "$file1" == "$file2" ]]; then
        printf "$file1\t$count1\t$count2\n" >> $output_dir/seq_count_summary.txt
    fi
    #Report file where no sequences were reoriented (i.e. the output was 0)
    grep -Fq $file1 tempdir2/seq_count_after.txt
    if [[ $? != 0 ]]; then
        printf "$file1\t$count1\t0\n" >> $output_dir/seq_count_summary.txt
    fi
done < tempdir2/seq_count.txt 

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
    rm *.$extension
fi

#Delete tempdir2
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir2" ]]; then
        rm -rf tempdir2
    fi
fi
}


##################################################################################
### Cleaning up and compiling final stats file UNIVERSAL fastx (but not for PE assembly and demux and GeneX) ###
##################################################################################
function clean_and_make_stats () {
countstart=$(date +%s)

#Delete empty output files
find $output_dir -empty -type f -delete
### Count reads before and after the process
touch tempdir2/seq_count_before.txt
seqkit stats --threads 6 -T *.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> tempdir2/seq_count_before.txt
touch tempdir2/seq_count_after.txt
outfile_check=$(ls $output_dir/*.$extension 2>/dev/null | wc -l)
if (( $outfile_check != 0 )); then
    seqkit stats --threads 6 -T $output_dir/*.$extension | awk -F'\t' 'BEGIN{OFS="\t";} NR!=1 {print $1,$4}' >> tempdir2/seq_count_after.txt
else
    printf '%s\n' "ERROR]: no output files generated ($output_dir). Check settings!" >&2
    end_process
fi


### Compile a track reads summary file (seq_count_summary.txt)
output_dir_for_sed=$(echo $output_dir | sed -e "s/\//\\\\\//g")
sed -e "s/$output_dir_for_sed\///" < tempdir2/seq_count_after.txt > tempdir2/seq_count_after.temp
printf "File\tReads_in\tReads_out\n" > $output_dir/seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    while read LINE2; do
        file2=$(echo $LINE2 | awk '{print $1}')
        count2=$(echo $LINE2 | awk '{print $2}')
        if [[ "$file1" == "$file2" ]]; then
            printf "$file1\t$count1\t$count2\n" >> $output_dir/seq_count_summary.txt
        fi
    done < tempdir2/seq_count_after.temp
    #Report file where no sequences were reoriented (i.e. the output was 0)
    grep -Fq $file1 tempdir2/seq_count_after.temp
    if [[ $? != 0 ]]; then
        printf "$file1\t$count1\t0\n" >> $output_dir/seq_count_summary.txt
    fi
done < tempdir2/seq_count_before.txt

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
    rm *.$extension
fi
#Remove mothur logfiles
mothur_logfiles=$(ls -1 *.logfile 2>/dev/null | wc -l)
if (( $mothur_logfiles != 0 )); then 
    rm mothur.*.logfile 
fi
#Delete tempdir
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir" ]]; then
        rm -rf tempdir
    fi
else 
    #compress files in /tempdir
    if [[ -d "tempdir" ]]; then
        pigz tempdir/*
    fi
fi
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir2" ]]; then
        rm -rf tempdir2
    fi
fi
}

##########################################################
### Cleaning up and compiling final stats file, when outputting multiple directories ###
##########################################################
function clean_and_make_stats_multidir () {
### Count reads before and after the process
mkdir -p tempdir2

if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    touch tempdir2/seq_count_before.txt
    for file in *.$extension; do
        size=$(echo $(cat $file | wc -l) / 4 | bc)
        printf "$file\t$size\n" >> tempdir2/seq_count_before.txt
    done
    touch tempdir2/seq_count_after.txt
    outfile_check=$(ls $output_dir/$subdir/*.$extension 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then 
        for file in $output_dir/$subdir/*.$outfile_addition.$extension; do
            size=$(echo $(cat $file | wc -l) / 4 | bc)
            printf "$file\t$size\n" >> tempdir2/seq_count_after.txt
        done
    fi
fi
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    touch tempdir2/seq_count_before.txt
    for file in *.$extension; do
        size=$(grep -c "^>" $file)
        printf "$file\t$size\n" >> tempdir2/seq_count_before.txt
    done
    touch tempdir2/seq_count_after.txt
    outfile_check=$(ls $output_dir/$subdir/*.$extension 2>/dev/null | wc -l)
    if (( $outfile_check != 0 )); then 
        for file in $output_dir/$subdir/*.$outfile_addition.$extension; do
            size=$(grep -c "^>" $file)
            printf "$file\t$size\n" >> tempdir2/seq_count_after.txt
        done
    fi
fi
### Compile a track reads summary file (seq_count_summary.txt)
output_dir_for_sed=$(echo $output_dir | sed -e "s/\//\\\\\//g")
sed -e "s/\.$outfile_addition//" < tempdir2/seq_count_after.txt | \
sed -e "s/^$output_dir_for_sed\///" | sed -e "s/^$subdir\///" > tempdir2/seq_count_after.temp
subdir=$(echo $subdir | sed -e "s/\\\\//g")
printf "File\tReads_in\tReads_out\n" > $output_dir/$subdir/seq_count_summary.txt
while read LINE; do
    file1=$(echo $LINE | awk '{print $1}')
    count1=$(echo $LINE | awk '{print $2}')
    while read LINE2; do
        file2=$(echo $LINE2 | awk '{print $1}')
        count2=$(echo $LINE2 | awk '{print $2}')
        if [[ "$file1" == "$file2" ]]; then
            printf "$file1\t$count1\t$count2\n" >> $output_dir/$subdir/seq_count_summary.txt
        fi
    done < tempdir2/seq_count_after.temp
    #Report file where no sequences were reoriented (i.e. the output was 0)
    grep -Fq $file1 tempdir2/seq_count_after.temp
    if [[ $? != 0 ]]; then
        printf "$file1\t$count1\t0\n" >> $output_dir/$subdir/seq_count_summary.txt
    fi
done < tempdir2/seq_count_before.txt
rm tempdir2/seq_count_after.txt tempdir2/seq_count_before.txt

#Note for counting seqs in FASTQ files
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    printf "\nPlease note that sequence count assumes that there are 4 lines per sequence in a FASTQ file (as this is mostly the case).
You may double-check the sequence count of one file using implemented 'QualityCheck' module in PipeCraft.\n" >> $output_dir/$subdir/seq_count_summary.txt
fi

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
    rm *.$extension
fi
#Remove mothur logfiles
mothur_logfiles=$(ls -1 *.logfile 2>/dev/null | wc -l)
if (( $mothur_logfiles != 0 )); then 
    rm mothur.*.logfile 
fi
#Delete tempdir
if [[ $debugger != "true" ]]; then
    if [[ -d "tempdir" ]]; then
        rm -rf tempdir
    fi
else 
    #compress files in /tempdir
    if [[ -d "tempdir" ]]; then
        pigz tempdir/*
    fi
fi
}


#################################################
### Paired-end data reorient reads based on FWD primers ###
#################################################
function PE_reorient_FWD () {
touch tempdir/R1.5_3.fastq
touch tempdir/R2.5_3.fastq
### Reorient based on FWD primer(s)
for primer in $(echo $fwd_tempprimer | sed "s/,/ /g"); do
    printf '%s\n' " searching FWD primer $primer"
    #convert IUPAC codes in fwd primer
    fwd_primer=$(convert_IUPAC $primer)
    #seach fwd primer in R1 and write to file
    fqgrep -m $mismatches -p $fwd_primer -e $inputR1.$extension >> tempdir/R1.5_3.fastq
    #seach fwd primer in R2 and write to file
    fqgrep -m $mismatches -p $fwd_primer -e $inputR2.$extension >> tempdir/R2.5_3.fastq
done
}

#################################################
### Paired-end data reorient reads based on REV primers ###
#################################################
function PE_reorient_REV () {
touch tempdir/R1.3_5.fastq
touch tempdir/R2.3_5.fastq
### Reorient based on REV primer(s)
for primer in $(echo $rev_tempprimer | sed "s/,/ /g"); do
    #convert IUPAC codes in rev primer
    rev_primer=$(convert_IUPAC $primer)
    printf '%s\n' " searching REV primer $primer"
    #search rev primer in R1 and write to file
    fqgrep -m $mismatches -p $rev_primer -e $inputR1.$extension >> tempdir/R1.3_5.fastq
    #search rev primer in R2 and write to file
    fqgrep -m $mismatches -p $rev_primer -e $inputR2.$extension >> tempdir/R2.3_5.fastq
done
}

#######################################################
### Single-end data reorient reads based on primers ###
#######################################################
#Fwd primers
function SE_reorient_FWD () {
touch tempdir/5_3.fastx
### Reorient based on FWD primer(s)
for primer in $(echo $fwd_tempprimer | sed "s/,/ /g"); do
    printf '%s\n' " searching FWD primer $primer"
    #convert IUPAC codes in fwd primer
    fwd_primer=$(convert_IUPAC $primer)
    #seach fwd primer and write to file
    if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
        fqgrep -m $mismatches -p $fwd_primer -e $input.$extension >> tempdir/5_3.fastx
    elif [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
        fqgrep -m $mismatches -p $fwd_primer -f -e $input.$extension >> tempdir/5_3.fastx
    else
        printf '%s\n' "ERROR]: $file formatting not supported!
Supported extensions: fastq, fq, fasta, fa, fas (and gz or zip compressed formats).
>Quitting" >&2
        end_process
    fi
done
}
#Rev primers
function SE_reorient_REV () {
touch tempdir/3_5.fastx
### Reorient based on REV primer(s)
for primer in $(echo $rev_tempprimer | sed "s/,/ /g"); do
    #convert IUPAC codes in rev primer
    rev_primer=$(convert_IUPAC $primer)
    printf '%s\n' " searching REV primer $primer"
    #search rev primer and write to file
    if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
        fqgrep -m $mismatches -p $rev_primer -e $input.$extension >> tempdir/3_5.fastx
    elif [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
        fqgrep -m $mismatches -p $rev_primer -f -e $input.$extension >> tempdir/3_5.fastx
    else
        printf '%s\n' "ERROR]: $file formatting not supported!
Supported extensions: fastq, fq, fasta, fa, fas (and gz or zip compressed formats).
>Quitting" >&2
        end_process
    fi
done
}

#############################
### Multiprimer artefacts ###
#############################
# R1 multi-primer artefacts search
function multiprimer_search_R1 () {
checkerror=$(seqkit rmdup --quiet -w 0 -n -D tempdir/duplicatesR1.temp tempdir/R1.5_3.fastq > tempdir/R1.5_3.fastq.temp 2>&1)
check_app_error
if [[ -s tempdir/duplicatesR1.temp ]]; then
    awk 'BEGIN{FS=","}{print $2}' tempdir/duplicatesR1.temp | sed -e 's/^ //' > tempdir/duplicatesR1.names
        #Remove duplicate seqs from fastq
    checkerror=$(seqkit grep --invert-match -n -w 0 -f tempdir/duplicatesR1.names tempdir/R1.5_3.fastq.temp \
    -o tempdir/$inputR1.$extension 2>&1)
    check_app_error
        #Get multi-primer artefacts
    checkerror=$(seqkit grep -w 0 -f tempdir/duplicatesR1.names tempdir/R1.5_3.fastq.temp -o tempdir/$inputR1.multiprimer.$extension 2>&1)
    check_app_error
    multiprimer_count=$(wc -l tempdir/duplicatesR1.names | awk '{print $1}')
    printf "   - found $multiprimer_count 'multi-primer' chimeric sequence(s) from $inputR1.$extension \n"
else
    mv tempdir/R1.5_3.fastq.temp tempdir/$inputR1.$extension
    printf "   - no 'multi-primer' chimeric sequences found from $inputR1.$extension \n"
fi
}

# R2 multi-primer artefacts search
function multiprimer_search_R2 () {
checkerror=$(seqkit rmdup --quiet -w 0 -n -D tempdir/duplicatesR2.temp tempdir/R2.3_5.fastq > tempdir/R2.3_5.fastq.temp 2>&1)
check_app_error
if [[ -s tempdir/duplicatesR2.temp ]]; then
    awk 'BEGIN{FS=","}{print $2}' tempdir/duplicatesR2.temp | sed -e 's/^ //' > tempdir/duplicatesR2.names
        #Remove duplicate seqs from fastq
    checkerror=$(seqkit grep --invert-match -n -w 0 -f tempdir/duplicatesR2.names tempdir/R2.3_5.fastq.temp \
    -o tempdir/$inputR2.$extension 2>&1)
    check_app_error
        #Get multi-primer artefacts
    checkerror=$(seqkit grep -w 0 -f tempdir/duplicatesR2.names tempdir/R2.3_5.fastq.temp -o tempdir/$inputR2.multiprimer.$extension 2>&1)
    check_app_error

    multiprimer_count=$(wc -l tempdir/duplicatesR2.names | awk '{print $1}')
    printf "   - found $multiprimer_count 'multi-primer' chimeric sequence(s) from $inputR2.$extension \n"
else
    mv tempdir/R2.3_5.fastq.temp tempdir/$inputR2.$extension
    printf "   - no 'multi-primer' chimeric sequences found from $inputR2.$extension \n"
fi
}

# Single-end data multi-primer artefacts search
function multiprimer_search_SE () {
checkerror=$(seqkit rmdup --quiet -w 0 -n -D tempdir/duplicates.temp tempdir/5_3.fastx > tempdir/5_3.fastx.temp 2>&1)
check_app_error
if [[ -s tempdir/duplicates.temp ]]; then
    awk 'BEGIN{FS=","}{print $2}' tempdir/duplicates.temp | sed -e 's/^ //' > tempdir/duplicates.names
        #Remove duplicate seqs from fastx
    checkerror=$(seqkit grep --invert-match -n -w 0 -f tempdir/duplicates.names tempdir/5_3.fastx.temp \
    -o tempdir/$input.$extension 2>&1)
    check_app_error

        #Get multi-primer artefacts
    checkerror=$(seqkit grep -w 0  -f tempdir/duplicates.names tempdir/5_3.fastx.temp \
    -o tempdir/$input.multiprimer.$extension 2>&1)
    check_app_error

    multiprimer_count=$(wc -l tempdir/duplicates.names | awk '{print $1}')
    printf "   - found $multiprimer_count 'multi-primer' chimeric sequence(s) from $input.$extension \n"
else
    mv tempdir/5_3.fastx.temp tempdir/$input.$extension
    printf "   - no 'multi-primer' chimeric sequences found from $input.$extension \n"
fi
}


###########################
### Check barcodes file ###
###########################
### Check barcodes file
function check_indexes_file () {
printf "\nValidating indexes file ...\n"
    #is fasta format?
cat $indexes_file | seqkit seq -v -w 0 > tempdir2/ValidatedBarcodesFileForDemux.fasta.temp
if [[ "$?" != "0" ]]; then
    printf '%s\n' "ERROR]: 'indexes file' not in correct fasta format. 
Please check the indexes file and format according to the 'indexes_file_example.txt'
>Quitting" >&2
    end_process
fi
    #does not contain duplicate values?
sort tempdir2/ValidatedBarcodesFileForDemux.fasta.temp | \
uniq --count --repeated | \
grep . && printf '%s\n' "ERROR]: indexes or samples names (above) are not unique in the indexes file. 
Please check the indexes file and include only unique index combinations and sample names.
>Quitting" >&2 && end_process \
|| :
    # report if dual indexes or single indexes
if grep -q "\..." tempdir2/ValidatedBarcodesFileForDemux.fasta.temp; then
    echo "ok; dual indexes"
    tag=$"dual"
    export tag
else
    echo "ok; single indexes"
    tag=$"single"
    export tag
fi
}


#######################################################
### Function to count features and sequences in OTU/ASV table ###
#######################################################
 # Handles "Sequence" column in the table
    # output variables:
    # feature_count - numbers of ASVs/OTUs in the table
    # nSeqs - numbers of sequences in the table
    # nSample - numbers of samples in the table
function count_features() {
    local table=$1
    
    # count ASVs/OTUs (rows minus header)
    feature_count=$(awk 'NR>1' "$table" | wc -l)
    
    # count sequences and samples
    while IFS='=' read -r key value; do
        case $key in
            "samples") nSample=$value ;;
            "sequences") nSeqs=$value ;;
        esac
    done < <(awk -F'\t' '
    NR==1 {
        seq_col = -1
        for(i=1; i<=NF; i++) {
            if($i == "Sequence") {
                seq_col = i
                break
            }
        }
        nSample = NF - 1 - (seq_col > 0 ? 1 : 0)
        print "samples=" nSample
    }
    NR>1 {
        sum = 0
        for(i=2; i<=NF; i++) {
            if(i != seq_col) {
                sum += $i
            }
        }
        total_seqs += sum
    }
    END {
        print "sequences=" total_seqs
    }' "$table")
    
    export feature_count
    export nSeqs
    export nSample
}
  

################################################
### generate README.txt for table filtering (curate_table_wf.sh) ###
################################################
function readme_table_filtering() {
    local output_dir=$1
    local runtime=$2

    # Start README with basic info
    cat << EOF > "$output_dir/README.txt"
# Feature table curation.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input parameters:
---------------
Tag-jumps filtering: ${filter_tag_jumps}
- f-value: ${f_value}
- p-value: ${p_value}
Length filtering:
- min length: ${min_length_num} (0 means OFF)
- max length: ${max_length_num} (0 means OFF)
Collapse identical ASVs/OTUs: ${collapseNoMismatch}

Output files:
------------
EOF

    # Add file descriptions based on what was generated
    if [[ $filter_tag_jumps == "true" ]]; then
        count_features "$output_dir/${feature_table_base_name%%.txt}_TagJumpFilt.txt"

        cat << EOF >> "$output_dir/README.txt"
## Tag-jumps filtering output (applied first):
   - *_TagJumpFilt.txt = Feature table after tag-jumps filtering
   - $fasta_base_name = Representative sequences after tag-jumps filtering (number of seqs = $ASVs_count, here always same as input)
   - TagJump_plot.pdf = Visualization of tag-jumps based on parameters
   - TagJump_stats.txt = Statistics about tag-jumps filtering including:
                             * Total reads
                             * Number of tag-jump events
                             * Tag-jump reads
                             * Read percent removed
   - tag-jumps_filt.log = R log file for tag-jumps filtering

    Number of Features                       = $feature_count
    Number of sequences in the Feature table = $nSeqs
    Number of samples in the Feature table   = $nSample

EOF
    fi

    if [[ -n "$feature_table_file2" ]]; then
        echo "   [clustered zOTUs (OTU_table) was also tag-jumps filtered, but stats not given in this file]" >> "$output_dir/README.txt"
        echo "" >> "$output_dir/README.txt"
    fi

    if [[ $collapseNoMismatch == "true" ]]; then
        count_features "$output_dir/${feature_table_base_name%%.txt}_collapsed.txt"
        cat << EOF >> "$output_dir/README.txt"
## Collapsing identical ASVs/OTUs output (applied after tag-jumps filtering, if ON):
[before collapsing, length filter is also applied (if ON)]
$ASVs_collapsed_result

EOF
     if [[ $feature_count > 0 ]]; then
     cat << EOF >> "$output_dir/README.txt"
        Number of Features                       = $feature_count
        Number of sequences in the Feature table = $nSeqs
        Number of samples in the Feature table   = $nSample
EOF
    fi 
    
    elif [[ $collapseNoMismatch == "false" ]] && \
         [[ ($min_length_num != "0" && -n $min_length_num) || \
            ($max_length_num != "0" && -n $max_length_num) ]]; then
        
        cat << EOF >> "$output_dir/README.txt"
## Length filtering output (applied after tag-jumps filtering, if ON):
$ASVs_lenFilt_result

EOF

    if [[ -n "$feature_table_file2" ]]; then
        echo "  [clustered zOTUs (OTU_table) was also length filtered (after tag-jumps filtering), but stats not given in this file]" >> "$output_dir/README.txt"
    fi


    fi

    # Add citations
    cat << EOF >> "$output_dir/README.txt"

##############################################
### Third-party applications for this process:
##############################################
# vsearch (version $vsearch_version)
   Citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584 https://doi.org/10.7717/peerj.2584
   https://github.com/torognes/vsearch
# seqkit (version $seqkit_version)
   Citation: Shen W, Le S, Li Y, Hu F (2016) 
   SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. 
   PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
# UNCROSS2
   Citation: R.C. Edgar (2018), UNCROSS2: identification of cross-talk in 16S rRNA OTU tables, https://doi.org/10.1101/400762
# R (version $R_version)
   Citation: R Core Team (2023). R: A language and environment for statistical computing. 
   R Foundation for Statistical Computing, Vienna, Austria.
   https://www.R-project.org/
##############################################
EOF
}

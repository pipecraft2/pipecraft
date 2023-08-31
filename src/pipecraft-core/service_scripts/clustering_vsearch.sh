#!/bin/bash

# Sequence clustering with vsearch
#Input = single-end fasta/fastq files.
#Output = FASTA formated representative OTU sequences and OTU_table.txt.

##########################################################
###Third-party applications:
#vsearch v2.23.0
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422
    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
    #Copyright (C) 2007-2021 Ole Tange, http://ole.tange.dk and Free Software Foundation, Inc.
    #Distributed under the License GPLv3+
#pigz v2.4
##########################################################
#load variables
id=$"--id ${similarity_threshold}"          # positive float (0-1)
otutype=$"--${OTU_type}"                    # list: --centroids, --consout
strands=$"--strand ${strands}"              # list: --strand both, --strand plus
remove_singletons=$"${remove_singletons}"   # true/false

#additional options
seqsort=$"${sequence_sorting}"           # list: --cluster_size or --cluster_fast, --cluster_smallmem
simtype=$"--iddef ${similarity_type}"    # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
centroid=$centroid_type                  # list: similarity, abundance
maxaccepts=$"--maxaccepts ${maxaccepts}" # pos integer
mask=$"--qmask ${mask}"                  # list: --qmask dust, --qmask none
cores=$"--threads ${cores}"              # pos integer
###############################
# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/clustering_out"

#additional options, if selection != undefined/false
if [[ $seqsort == "size" ]]; then
    seqsort=$"--cluster_size"
elif [[ $seqsort == "length" ]]; then
    seqsort=$"--cluster_fast"
elif [[ $seqsort == "none" ]]; then
    seqsort=$"--cluster_smallmem --usersort"
fi 
if [[ $centroid == "similarity" ]]; then
    centroid_in=$"" 
else
    centroid_in=$"--sizeorder"
fi

#############################
### Start of the workflow ###
#############################
start=$(date +%s)

### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env

### Pre-process samples
printf "Checking files ... \n"
for file in *.$fileFormat; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_gz_zip_SE
    ### Check input formats
    check_extension_fastx
done

#If input is FASTQ then convert to FASTA
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    for file in *.$extension; do
        samp_name=$(basename $file | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
        checkerror=$(seqkit fq2fa -t dna --line-width 0 $file -o $samp_name.fasta 2>&1)
        check_app_error
    done

    was_fastq=$"true"
    export was_fastq
    extension=$"fasta"
    export extension
fi

#tempdir
if [[ -d tempdir ]]; then
    rm -rf tempdir
fi
mkdir -p tempdir

### Rename sequences to sha1
    # and dereplication of individual samples, add sample ID to the header
derep_rename () {
  samp_name=$(basename $1 | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
  vsearch \
    --derep_fulllength "$1" \
    --relabel_sha1 \
    --output - \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout \
  | sed 's/>.*/&;sample='"$samp_name"'/' > tempdir/"$samp_name".fasta
}
export -f derep_rename
printf "Dereplication of individual samples ... \n"
find . -maxdepth 1 -name "*.$extension" | parallel -j 1 "derep_rename {}"

### Global dereplication
printf "Dereplicating globally ... \n"
find tempdir -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
| vsearch \
--derep_fulllength - \
--output $output_dir/Glob_derep.fasta \
--uc tempdir/Glob_derep.uc \
--fasta_width 0 \
--threads 1 \
--sizein --sizeout

### Clustering
printf "Clustering ... \n"
printf "\n vsearch $seqsort \
$output_dir/Glob_derep.fasta \
$id \
$simtype \
$strands \
$mask \
$centroid_in \
$maxaccepts \
$cores \
$otutype $output_dir/OTUs.temp.fasta \
--uc $output_dir/OTUs.uc \
--fasta_width 0 \
--sizein "

checkerror=$(vsearch $seqsort \
$output_dir/Glob_derep.fasta \
$id \
$simtype \
$strands \
$mask \
$centroid_in \
$maxaccepts \
$cores \
$otutype $output_dir/OTUs.temp.fasta \
--uc $output_dir/OTUs.uc \
--fasta_width 0 \
--sizein --sizeout 2>&1)
check_app_error

### Cat dereplicated individual samples for making an OTU table
cat tempdir/*.fasta > tempdir/Dereplicated_samples.fasta

## Prepare table with sequence abundance per sample
seqkit seq --name tempdir/Dereplicated_samples.fasta \
  | awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
  | sed 's/size=//; s/sample=//' \
  > tempdir/ASV_table_long.txt

### OTU table creation
printf "Making OTU table ... \n"
Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
  --derepuc      tempdir/Glob_derep.uc \
  --uc           "$output_dir"/OTUs.uc \
  --asv          tempdir/ASV_table_long.txt \
  --rmsingletons $remove_singletons \
  --output       "$output_dir"/OTU_table.txt 2>&1)
echo $Rlog > tempdir/OTU_table_creation.log 
wait

### Discard singleton OTUs
if [[ $remove_singletons == "true"  ]]; then
    printf "Discarding singletons ... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --sizeout --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error
    # remove ";sample=.*;" from OTU.fasta file.
    sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
    # removing ";size=" because OTU table does not have "size" annotations; so the files would fit to LULU
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta 
    mv $output_dir/OTUs.temp.fasta tempdir/
else
    sed -e 's/;sample=.*;/;/' $output_dir/OTUs.temp.fasta > $output_dir/OTUs.fasta
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta
    mv $output_dir/OTUs.temp.fasta tempdir/
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

mv $output_dir/Glob_derep.fasta tempdir/Glob_derep.fasta

#If input = FASTQ, then mkdir for converted fasta files
if [[ $was_fastq == "true" ]]; then
    mkdir -p $output_dir/clustering_input_to_FASTA
    mv *.fasta $output_dir/clustering_input_to_FASTA
fi

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
  delete_uncompressed=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $(NF-1)}';))
  rm *.$delete_uncompressed
fi

#Delete tempdirs
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir ]]; then
        rm -rf tempdir
    fi
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    if [[ -f $output_dir/OTU_table_creation.log ]]; then
        rm -f $output_dir/OTU_table_creation.log
    fi
    else 
    #compress files in /tempdir
    if [[ -d tempdir ]]; then
        pigz tempdir/*
    fi
fi

#Make README.txt file
OTU_count=$(grep -c "^>" $output_dir/OTUs.fasta)
nSeqs=$(awk 'BEGIN{FS=OFS="\t"}NR>1{for(i=2;i<=NF;i++) t+=$i; print t; t=0}' $output_dir/OTU_table.txt | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
nCols=$(awk -F'\t' '{print NF; exit}' $output_dir/OTU_table.txt)
nSample=$(awk -v NUM=$nCols 'BEGIN {print (NUM-1)}') # -1 cuz 1st column is OTU_ID

end=$(date +%s)
runtime=$((end-start))

printf "# Reads were clustered to OTUs using vsearch (see 'Core command' below for the used settings).

Number of OTUs                       = $OTU_count
Number of sequences in the OTU table = $nSeqs
Number of samples in the OTU table   = $nSample

Files in 'clustering_out' directory:
# OTUs.fasta    = FASTA formated representative OTU sequences. OTU headers are renamed according to sha1 algorithm in vsearch.
# OTU_table.txt = OTU distribution table per sample (tab delimited file). OTU headers are renamed according to sha1 algorithm in vsearch.
# OTUs.uc       = uclust-like formatted clustering results for OTUs.

Core command -> 
clustering: vsearch $seqsort dereplicated_sequences.fasta $id $simtype $strands $mask $centroid_in $maxaccepts $cores $otutype OTUs.fasta " > $output_dir/README.txt

## if input was fastq
if [[ $was_fastq == "true" ]]; then
  printf "\n\nInput was fastq; converted those to fasta before clustering. 
  Converted fasta files in directory 'clustering_input_to_FASTA' \n" >> $output_dir/README.txt
fi

printf "\nTotal run time was $runtime sec.\n\n
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.23.0 for clustering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422 for job parallelisation 
    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
##########################################################" >> $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

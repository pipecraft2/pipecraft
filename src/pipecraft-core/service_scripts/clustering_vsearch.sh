#!/bin/bash

# Sequence clustering with vsearch
#Input = single-end fasta/fastq files.
#Output = FASTA formated representative OTU sequences and OTU_table.txt.

##########################################################
###Third-party applications:
#vsearch v2.22.1
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
echo "WORKDIR: $workingDir"

#load variables
extension=$fileFormat
#mandatory options
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
if [[ $remove_singletons == "true"  ]]; then
    remove_singletons=$"TRUE"
fi
if [[ $remove_singletons == "false"  ]]; then
    remove_singletons=$"FALSE"
fi

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check_clustering
### Prepare working env and check paired-end data
prepare_SE_env

### Pre-process samples
printf "Checking files ... \n"
for file in *.$extension; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$extension//")
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        
    check_gz_zip_SE
    ### Check input formats (only fasta, fa, fas supported)
    check_extension_fasta
done

#tempdir
if [ -d tempdir ]; then
    rm -rf tempdir
fi
mkdir -p tempdir

### Rename sequences to md5
    # and dereplication of individual samples, add sample ID to the header
derep_rename () {
  samp_name=$(basename $1 | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
  vsearch \
    --derep_fulllength "$1" \
    --relabel_md5 \
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
--sizein --sizeout"

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
echo $Rlog > $output_dir/R_run.log 
wait

### Discard singleton OTUs
if [[ $remove_singletons == "TRUE"  ]]; then
    printf "Discarding singletons ... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --sizeout --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error

    sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
else
    sed -e 's/;sample=.*;/;/' $output_dir/OTUs.temp.fasta > $output_dir/OTUs.fasta
    rm $output_dir/OTUs.temp.fasta
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
    rm *.$extension
fi

#Delete tempdirs
if [[ -d tempdir ]]; then
    rm -rf tempdir
fi
if [[ -d tempdir2 ]]; then
    rm -rf tempdir2
fi
rm $output_dir/Glob_derep.fasta
#rm
if [[ -f $output_dir/R_run.log ]]; then
    rm -f $output_dir/R_run.log
fi

size=$(grep -c "^>" $output_dir/OTUs.fasta)
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "Clustering formed $size OTUs.

Files in 'clustering_out' directory:
# OTUs.fasta    = FASTA formated representative OTU sequences. OTU headers are renamed according to MD5 algorithm in vsearch.
# OTU_table.txt = OTU distribution table per sample (tab delimited file). OTU headers are renamed according to MD5 algorithm in vsearch.
# OTUs.uc       = uclust-like formatted clustering results for OTUs.

Core commands -> 
clustering: vsearch $seqsort dereplicated_sequences.fasta $id $simtype $strands $mask $centroid_in $maxaccepts $cores $otutype OTUs.fasta --fasta_width 0 --sizein --sizeout

Total run time was $runtime sec.\n\n
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.18.0 for clustering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422 for job parallelisation 
    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
##########################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Check README.txt files in output directory for further information about the process.\n"


printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "dataFormat=$dataFormat"
echo "readType=single_end"

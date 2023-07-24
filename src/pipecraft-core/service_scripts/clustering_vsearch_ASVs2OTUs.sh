#!/bin/bash

# ASVs to OTUs with vsearch
#Input = ASV fasta and ASV table
    # size annotation is important for ASVs but if not provided, then getting this info from the ASV table.
    # ASV table format = ASVs in rows, samples in columns (2nd col may be Sequence col (will be removed on the background))
#Output = FASTA formated representative OTU sequences and OTU table

##########################################################
###Third-party applications:
#vsearch v2.22.1
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
##########################################################
echo "WORKDIR: $workingDir"

#load data
ASV_fasta=${ASV_fasta}
#If specified, get OTUs subset
if [[ $ASV_fasta == "undefined" ]]; then
    :
else
    regex='[^/]*$'
    ASV_fasta_temp=$(echo $ASV_fasta | grep -oP "$regex")
    ASV_fasta=$(printf "/extraFiles/$ASV_fasta_temp")
    printf "\n ASV_fasta = $ASV_fasta \n"
fi
extension=$(basename $ASV_fasta | awk 'BEGIN{FS="."}{print $NF}')

ASV_table=${ASV_table}
#If specified, get OTUs subset
if [[ $ASV_fasta == "undefined" ]]; then
    :
else
    regex='[^/]*$'
    ASV_table_temp=$(echo $ASV_table | grep -oP "$regex")
    ASV_table=$(printf "/extraFiles2/$ASV_table_temp")
    printf "\n ASV_table = $ASV_table \n"
fi

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

### Check input formats (only fasta, fa, fas supported)
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported!
Supported extensions: fasta, fas, fa.
>Quitting" >&2
    end_process
fi

#tempdir
if [ -d tempdir ]; then
    rm -rf tempdir
fi
mkdir -p tempdir


### Check ASV table
    #and remove ;size= from the table
if [[ $(grep -c ";size" <(head $ASV_table)) > 0 ]]; then
    printf "\nRemoving size annotations from the ASV table\n"
    sed -e 's/;size=[^\t]*//' $ASV_table > tempdir/ASV_table.temp
else
    cp $ASV_table tempdir/ASV_table.temp
fi 
    #and drop 2nd col if it's a sequence
a=$(seqkit seq --validate-seq <(head -2 <(awk 'BEGIN{FS="\t"}{print $2}' tempdir/ASV_table.temp) | sed "1 s/^/>/"))
if [[ $? == "0" ]]; then
    printf "\nRemoving sequence column from the ASV table\n"
    awk 'BEGIN {FS=OFS="\t"} {$2=""; gsub("\t+", "\t", $0)}1' tempdir/ASV_table.temp > tempdir/ASV_table.txt && rm tempdir/ASV_table.temp
fi 

### Check size annotations in ASV_fasta  ----- THIS TAKSE FOREVER HERE!!! SHITFUCK - try this in R
if [[ $(grep -c ";size" <(head $ASV_fasta)) > 0 ]]; then
    echo "fasta size annotations - OK"
    cp $ASV_fasta tempdir/ASV_fasta.fasta
else 
    echo "no fasta size annotations - getting those from the ASV table"
    awk 'NR>1{for(i=1;i<=NF;i++) t+=$i; print $1";size="t; t=0}' tempdir/ASV_table.txt > tempdir/ASV_fasta.size

    cp $ASV_fasta tempdir/ASV_fasta.fasta
    time while read NAME; do
        id=$(echo $NAME | sed -e "s/;size=.*//")
        size=$(echo $NAME | sed -e "s/.*;size=//")
        
        sed -i "s/\<$id\>/$id;size=$size/" tempdir/ASV_fasta.fasta
    done < tempdir/ASV_fasta.size
fi






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
    rm *.$newextension
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
echo "fileFormat=$newextension"
echo "dataFormat=$dataFormat"
echo "readType=single_end"

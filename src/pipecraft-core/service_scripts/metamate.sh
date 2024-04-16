#!/bin/bash

# metaMATE (https://github.com/tjcreedy/metamate) software for 
#  removel of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data

# Input = fasta file of the OTUs/ASVs and the OTU/ASV table. 

##########################################################
###Third-party applications:
# metaMATE
    #citation: Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A.J., Vogler, A.P. and Emerson, B.C. (2021), Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data. Mol Ecol Resour, 21: 1772-1787. https://doi.org/10.1111/1755-0998.13337
    #Distributed under the GNU General Public License
    #https://github.com/tjcreedy/metamate
##########################################################

#load variables
find_or_dump=${find_or_dump}       # find, dump or find_and_dump
specifications=${specifications}   # file
reference_seqs=${reference_seqs}   # file
table=${table}                     # file 
rep_seqs=${rep_seqs}               # file
genetic_code=${genetic_code}       # integer
length=${length}                   # integer
find_results=${find_results}       # file
result_index=${result_index}       # integer
NA_abund_thresh=${NA_abund_thresh} # float
base_variation=${base_variation}   # integer
cores=${cores}                     # integer

# modify the path to the input files
regex='[^/]*$'
specifications=$(echo $specifications | grep -oP "$regex")
specifications=$(basename $specifications) #basename, needed for macOS
specifications=$(printf "/extraFiles2/$specifications")

reference_seqs=$(echo $reference_seqs | grep -oP "$regex")
reference_seqs=$(basename $reference_seqs)
reference_seqs=$(printf "/extraFiles3/$reference_seqs")

table=$(echo $table | grep -oP "$regex")
table=$(basename $table)
table=$(printf "/extraFiles4/$table")

rep_seqs=$(echo $rep_seqs | grep -oP "$regex")
rep_seqs=$(basename $rep_seqs)
rep_seqs=$(printf "/extraFiles5/$rep_seqs")

find_results=$(echo $find_results | grep -oP "$regex")
find_results=$(basename $find_results)
find_results=$(printf "/extraFiles8/$find_results")

printf "\n specifications file = $specifications\n"
printf "reference seqs file = $reference_seqs\n"
printf "table file = $table\n"
printf "rep seqs file = $rep_seqs\n"
printf "find results file = $find_results\n"
###############################
# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/metamate_out"

# activate metamate conda env in the container 
eval "$(conda shell.bash hook)"
conda activate metamate


# if perfoming clade binning, then ERROR when processing more than 65,536 ASVs
ASVcount=$(grep -c "^>" $rep_seqs)
if (( $ASVcount > 65536 )); then
    
fi

# FIND
# python3 metamate.py find \
#     -A $rep_seqs \
#     -M $table \
#     -S $specifications \
#     -R $reference_seqs \
#     --expectedlength 418 \
#     --basevariation 0 \
#     --table $genetic_code \
#     -t $cores \
#     -o $output_dir \
#     --overwrite
    
# DUMP 
# python3 metamate.py dump \
#     -A $rep_seqs \
#     -C metamate_find/resultcache \
#     -o metamate_dump/output \
#     --overwrite \
#     -i $result_index


#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=paired_end"

#!/bin/bash

# metaMATE v0.4.0 (https://github.com/tjcreedy/metamate) software for 
 #  removel of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data

# Input = fasta file of the OTUs/ASVs and the OTU/ASV table.

##########################################################
###Third-party applications:
# metaMATE v0.4.0
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


# if perfoming clade binning, then WARNING when processing more than 65,536 ASVs
ASVcount=$(grep -c "^>" $rep_seqs)
if (( $ASVcount > 65536 )); then
    printf '%s\n' "WARNING]: clade binning NOT performed, because the input ASVs limit is 65,536 for that.
Current input has $ASVcount ASVs."
fi

# /Files/metamate/metamate/metamate.py

#############################
### Start of the workflow ###
#############################
start=$(date +%s)

# FIND
python3 /metamate/metamate/metamate.py find \
    -A $rep_seqs \
    -M $table \
    -S $specifications \
    -R $reference_seqs \
    --expectedlength $length \
    --basesvariation $base_variation \
    --table $genetic_code \
    -t $cores \
    -o $output_dir \
    --overwrite

# DUMP 
# python3 metamate.py dump \
#     -A $rep_seqs \
#     -C metamate_find/resultcache \
#     -o metamate_dump/output \
#     --overwrite \
#     -i $result_index


#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    rm $output_dir/${rep_seqs%.*}_ASVcounts.csv # rm metaMATE formatted ASV table
fi

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
#count ASVs
#ASV_count=$(grep -c "^>" $output_dir2/ASVs.fasta)
#nSeqs=$(awk 'BEGIN{FS=OFS="\t"}NR>1{for(i=2;i<=NF;i++) t+=$i; print t; t=0}' $output_dir2/ASVs_table.txt  | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
#nCols=$(awk -F'\t' '{print NF; exit}' $output_dir2/ASVs_table.txt)
#nSample=$(awk -v NUM=$nCols 'BEGIN {print (NUM-2)}') # -2 cuz 1st column is ASV_ID and 2nd is Sequence

if (( $ASVcount > 65536 )); then
    warn=$(echo "WARNING]: clade binning NOT performed, because the input ASVs limit is 65,536 for that. Current input has $ASVcount ASVs.")
fi

if [[ $find_or_dump == "find" ]] || [[ $find_or_dump == "find_and_dump" ]]; then

    printf "# 'find' function of metaMATE to detect putative NUMT and other erroneous sequences.

$warn

Files in 'metamate_out' directory:
# ${rep_seqs%.*}_aligned.fasta = aligned $rep_seqs
# ${rep_seqs%.*}_UPGMA.tre     = newick-format UPGMA tree file (if clade binning was specified)
# results.csv                  = metaMATE find results file, which synthesises all of the results from applying all combinations of the specified terms and thresholds to the input ASVs, given the control groups of authentic- and non-authentic-ASVs.
        see this metaMATE documentation link for more info about the results file: https://github.com/tjcreedy/metamate?tab=readme-ov-file#results-find-only 
# resultcache                  = needed for dump. This is a compressed text file containing information on the ASVs rejected or retained for each of the supplied specification terms and threshold sets of a find run.
# ${rep_seqs%.*}_control.txt   = a two-column tab-separated table recording all ASVs determined to be validated-authentic or validated-non-authentic. 
*_clades.csv                   = a two-column comma-separated table recording the clade grouping for each input ASV.
# -> more info: https://github.com/tjcreedy/metamate?tab=readme-ov-file#outputs 

# Check out analyse_results_draft.R provided by the metaMATE developers for generating plots from a metaMATE find. 
  https://github.com/tjcreedy/metamate/blob/main/analyse_results_draft.R

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#metaMATE v0.4.0
    #citation: Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A.J., Vogler, A.P. and Emerson, B.C. (2021), Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data. Mol Ecol Resour, 21: 1772-1787. https://doi.org/10.1111/1755-0998.13337
    #https://github.com/tjcreedy/metamate
########################################################" > $output_dir/README.txt
fi 


if [[ $find_or_dump == "dump" ]] || [[ $find_or_dump == "find_and_dump" ]]; then
    printf "# 'dump' function of metaMATE to discard putative NUMT and other erroneous sequences."
fi

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"
    
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
find_or_dump=${find_or_dump}   # find, dump or find_and_dump
specifications=${specifications} # file
reference_seqs=${reference_seqs} # file
table=${table} # file 
rep_seqs=${rep_seqs} # file
genetic_code=${genetic_code} # integer
length=${length} # integer

find_results=${find_results} # file
result_index=${result_index} # integer
NA_abund_thresh=${NA_abund_thresh} # float
base_variation=${base_variation} # integer


echo $find_or_dump
echo $specifications
echo $rep_seqs
echo $table

regex='[^/]*$'
specifications=$(echo $specifications | grep -oP "$regex")
specifications=$(basename $specifications) #basename, needed for macOS
specifications=$(printf "/extraFiles/$specifications")

printf "\n new specifications = $specifications\n"

# FIND
# python3 metamate.py find \
#     -A $rep_seqs \
#     -M $table \
#     -S $specifications \
#     -R $reference_seqs \
#     --expectedlength 418 \
#     --basevariation 0 \
#     --table 5 \
#     -t 16 \
#     -o metamate_find \
#     --overwrite
    
# DUMP 
# python3 metamate.py dump \
#     -A $rep_seqs \
#     -M $table \
#     -o metamate_dump/output \
#     --overwrite \
#     -i $result_index


#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=paired_end"

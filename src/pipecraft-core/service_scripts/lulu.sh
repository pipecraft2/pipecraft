#!/bin/bash

# Post-clustering of OTU/ASV table with LULU 
# Input = OTU table (2nc col may be 'Sequence') and correspoding fasta file

################################################
###Third-party applications:
#lulu v0.1.0
    #citation: Froslev, T.G., Kjoller, R., Bruun, H.H. et al. (2017) Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nat Commun 8, 1188.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/tobiasgf/lulu
#BLAST 2.12.0+
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. 
#vsearch v2.23.0
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
#pigz
##############################################

# Checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
blast_version=$(blastn -version 2>&1 | head -n 1 | awk '{print $2}')
printf "# Checking tool versions ...\n"
printf "# vsearch (version $vsearch_version)\n"
printf "# BLAST (version $blast_version)\n"

#env variables
workingDir=${workingDir}
extension=$fileFormat && export fileFormat 
#load variables
match_list_soft=${match_list_soft}
vsearch_similarity_type=${vsearch_similarity_type}
perc_identity=${perc_identity}
coverage_perc=${coverage_perc}
strands=${strands}
cores=${cores}
#(variables used in lulu.R)
min_ratio_type=${min_ratio_type}
min_ratio=${min_ratio}
min_match=${min_match}
min_rel_cooccurence=${min_rel_cooccurence}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/lulu_out"
if [ -d $output_dir ]; then
    rm -rf $output_dir
fi
mkdir -p $output_dir

#get specified input OTU table file
regex='[^/]*$'
otu_table_temp=$(echo $table | grep -oP "$regex")
otu_table=$(printf "/extraFiles/$otu_table_temp")
export otu_table # for lulu.R
printf "\n input table = $otu_table \n"

#get specified input fasta file
regex='[^/]*$'
input_fasta_temp=$(echo $fasta_file | grep -oP "$regex")
input_fasta=$(printf "/extraFiles2/$input_fasta_temp")
export input_fasta # for lulu.R
printf "\n input fasta = $input_fasta \n"

#############################
### Start of the workflow ###
#############################
#start time
start_time=$(date)
start=$(date +%s)

### Check if input table col.names contain #
firstline=$(awk 'NR==1 {print; exit}' $otu_table)
if echo "$firstline" | grep -q "#"; then
    printf '%s\n' "ERROR]: first line of input table has # sign.
Please remove and try again.
>Quitting" >&2
    end_process
fi

### Generate match list for LULU
if [[ $match_list_soft == "BLAST" ]]; then
    printf "\n#Making blast database from the input fasta \n"
    checkerror=$(makeblastdb -in $input_fasta -parse_seqids -dbtype nucl 2>&1)
    check_app_error

    printf "# Generating match list for lulu using BLASTn \n"
    checkerror=$(blastn -db $input_fasta \
            -outfmt '6 qseqid sseqid pident' \
            -out $output_dir/match_list.lulu \
            -qcov_hsp_perc $coverage_perc \
            -perc_identity $perc_identity \
            -query $input_fasta \
            -num_threads $cores 2>&1)
            check_app_error
fi

if [[ $match_list_soft == "vsearch" ]]; then
    printf "# Generating match list for lulu using vsearch \n"
    #convert perc_identity and coverage_perc for vsearch
    vsearch_perc_identity=$(awk "BEGIN {print $perc_identity/100}")
    vsearch_coverage_perc=$(awk "BEGIN {print $coverage_perc/100}")
    #run vsearch
    checkerror=$(vsearch --usearch_global $input_fasta \
            --db $input_fasta \
            --strand both --self \
            --id $vsearch_perc_identity \
            --iddef $vsearch_similarity_type \
            --userout $output_dir/match_list.lulu \
            --userfields query+target+id \
            --maxaccepts 0 \
            --query_cov $vsearch_coverage_perc \
            --threads $cores 2>&1)
            check_app_error
fi

#Check if match list was generated
if [[ -e "$output_dir/match_list.lulu" ]]; then
    printf '%s\n' "match list generation with $match_list_soft done"
else
    printf '%s\n' "ERROR]: match list generation with $match_list_soft for LULU failed" >&2
    end_process
fi

###Run LULU in R
printf "# Running lulu\n"
errormessage=$(Rscript /scripts/submodules/lulu.R 2>&1)
echo $errormessage > $output_dir/lulu_R_run.log 
wait
printf "\n LULU completed \n"

# format R log
sed -i 's/;; /\n/g' $output_dir/lulu_R_run.log
sed -i '/progress:/d' $output_dir/lulu_R_run.log



# get output table name 
lulu_table=$(basename $otu_table | sed 's/\.[^.]*$//')".lulu.txt"

# get output fasta name 
lulu_fasta=$(basename $input_fasta | sed 's/\.[^.]*$//')".lulu.fasta"

# Generate new OTUs.fasta file that excluded "discarded" OTUs by lulu
checkerror=$(seqkit grep --quiet \
                         --line-width 0 \
                         -f $output_dir/lulu_out_OTUids.txt \
                         $input_fasta > $output_dir/$lulu_fasta 2>&1)
check_app_error

########################################
### CLEAN UP AND COMPILE README FILE ###
########################################
for i in $input_fasta.n*; do
    [[ -f $i ]] || continue
    rm -f "$i"
done
if [[ $debugger != "true" ]]; then
    # if [[ -f $output_dir/lulu_R_run.log ]]; then
    #     rm -f $output_dir/lulu_R_run.log
    # fi
    if [[ -f $output_dir/lulu_out_OTUids.txt ]]; then
        rm -f $output_dir/lulu_out_OTUids.txt
    fi
fi

end=$(date +%s)
runtime=$((end-start))

###Make README.txt file
#match_list generation
if [[ $match_list_soft == "BLAST" ]]; then
    match_list_generation=$"makeblastdb -in $input_fasta -parse_seqids -dbtype nucl; blastn -db $input_fasta -outfmt '6 qseqid sseqid pident' -out match_list.lulu -qcov_hsp_perc $coverage_perc -perc_identity $perc_identity -query $input_fasta"
fi
if [[ $match_list_soft == "vsearch" ]]; then
    match_list_generation=$"vsearch --usearch_global $input_fasta --db $input_fasta --strand both --self --id $vsearch_perc_identity --iddef $vsearch_similarity_type --userout match_list.lulu --userfields query+target+id --maxaccepts 0 --query_cov $vsearch_coverage_perc"
fi
#count merged units and curated units
curated_units=$(grep -c "^>" $output_dir/$lulu_fasta)
merged_units=$(wc -l $output_dir/discarded_units.lulu | awk '{print $1}')
if (( $merged_units == 0 )); then
    info=$"No output table generated; OTU/ASV counts remained the same with the selected settings."
    rm $output_dir/discarded_units.lulu
    rm $output_dir/$lulu_fasta
    rm $output_dir/$lulu_table
else
    info=$"Output table consists of $curated_units Features (OTUs/ASVs)."
fi

printf "# Performed post-clustering with LULU (see 'Core commands' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Total of $merged_units Features (OTUs/ASVs) were merged.
$info

Files in 'lulu_out' directory:
------------------------------
# $lulu_table     = curated table in tab delimited txt format
# $lulu_fasta     = fasta file for the Features (OTUs/ASVs) in the curated table
# match_list.lulu        = match list file that was used by LULU to merge 'daughter molecular' units
# discarded_units.lulu   = Features (OTUs/ASVs) that were merged with other units based on specified thresholds)

Core commands -> 
match list for LULU (match_list.lulu): $match_list_generation
LULU in R: curated_result <- lulu::lulu(otutable_name, match_list.lulu, minimum_ratio_type = $min_ratio_type, minimum_ratio = $min_ratio, minimum_match = $min_match, minimum_relative_cooccurence = $min_rel_cooccurence)

Total run time was $runtime sec.

################################################
###Third-party applications:
#lulu v0.1.0
    #citation: Froslev, T.G., Kjoller, R., Bruun, H.H. et al. 2017. Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nat Commun 8, 1188.
#BLAST (version $blast_version) (if BLAST was used to make match_list.lulu)
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. 
#vsearch (version $vsearch_version) (if vsearch was used to make match_list.lulu)
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
##############################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

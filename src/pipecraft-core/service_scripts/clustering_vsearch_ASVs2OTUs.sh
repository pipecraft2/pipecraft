#!/bin/bash

# ASVs to OTUs with vsearch
#Input = ASV fasta and ASV table
    # ASVs.fasta must be without size annotations; but size annotation are important, current workflow gets those from the ASV table.
    # ASV table format = ASVs in rows, samples in columns; 2nd column must be SEQUENCE column (default PipeCraft2 ASVs workflow output)
#Output = FASTA formated representative OTU sequences and OTU table

##########################################################
###Third-party applications:
#vsearch v2.23.0
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
##########################################################
#load ASV_fasta
ASV_fasta=${ASV_fasta}
if [[ $ASV_fasta == "undefined" ]]; then
    printf '%s\n' "ERROR]: ASV_fasta not selected. >Quitting" >&2
    end_process
else
    regex='[^/]*$'
    ASV_fasta_temp=$(echo $ASV_fasta | grep -oP "$regex")
    ASV_fasta=$(printf "/extraFiles/$ASV_fasta_temp")
    printf "\n ASV_fasta = $ASV_fasta \n"
fi
extension=$(basename $ASV_fasta | awk 'BEGIN{FS="."}{print $NF}')

#load ASV_table
ASV_table=${ASV_table}
if [[ $ASV_table == "undefined" ]]; then
    printf '%s\n' "ERROR]: ASV_table not selected. >Quitting" >&2
    end_process
else
    regex='[^/]*$'
    ASV_table_temp=$(echo $ASV_table | grep -oP "$regex")
    ASV_table=$(printf "/extraFiles2/$ASV_table_temp")
    printf "\n ASV_table = $ASV_table \n"
fi

#settings
id=$"--id ${similarity_threshold}"          # positive float (0-1)
otutype=$"--${OTU_type}"                    # list: --centroids, --consout
strands=$"--strand ${strands}"              # list: --strand both, --strand plus
remove_singletons=$"${remove_singletons}"   # true/false
seqsort=$"${sequence_sorting}"              # list: --cluster_size or --cluster_fast, --cluster_smallmem
simtype=$"--iddef ${similarity_type}"       # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
centroid=$centroid_type                     # list: similarity, abundance
maxaccepts=$"--maxaccepts ${maxaccepts}"    # pos integer
mask=$"--qmask ${mask}"                     # list: --qmask dust, --qmask none
cores=$"--threads ${cores}"                 # pos integer
###############################
# Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/ASVs2OTUs_out"
export output_dir

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
echo "output_dir = $output_dir"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir
if [[ -d "tempdir" ]]; then
    rm -rf tempdir
fi
mkdir -p tempdir

### Check the ASV table
head -2 $ASV_table | awk 'NR>1{print $2}' | sed -e "s/^/>first_seq\n/" | seqkit seq --validate-seq -w 0
#2nd column must be the sequence column
if [[ $? != 0 ]]; then
    printf '%s\n' "ERROR]: 2nd column does not seem to be SEQUENCE in the $ASV_table.
    Herein, 2nd column must be the sequence column.
    >Quitting" >&2
    end_process
fi
# must not contain size annotations
if [[ $(grep -c ";size=" <(head $ASV_table)) > 0 ]]; then
    printf '%s\n' "ERROR]: found size annotations (;size=) in $ASV_table.
    Not supported; please remove.
    >Quitting" >&2
    end_process
else
    export ASV_table
fi

### Check the ASV fasta
    # ERROR if "";size=" in the fasta
if [[ $(grep -c ";size=" <(head $ASV_fasta)) > 0 ]]; then
    printf '%s\n' "ERROR]: found size annotations (;size=) in $ASV_fasta.
    Not supported; please remove.
    >Quitting" >&2
    end_process
else
    ### Get size annotations for ASVs for clustering (making ASVs.size.fasta from the ASV table where 2nd column is SEQUENCE).
        #Sum each ROW in the table and print ">" + 1st col + row sum as ";size=" + 2nd col (sequence)
        #NR>1 skips the first row, which are the col names. 
    fastasize=$(basename $ASV_fasta | awk 'BEGIN{FS=OFS="."}NF{NF -=1}1')
    awk 'NR>1{for(i=3;i<=NF;i++) t+=$i; print ">"$1";size="t"\n"$2; t=0}' $ASV_table > $output_dir/$fastasize.size.fasta
    #export for R
    ASV_fasta_size=$"$output_dir/$fastasize.size.fasta"
    export ASV_fasta_size
fi

### Clustering
checkerror=$(vsearch $seqsort \
                    $ASV_fasta_size \
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


### OTU table creation
printf "Making OTU table ... \n "
Rlog=$(Rscript /scripts/submodules/ASVs2OTUs.R 2>&1)
echo $Rlog > $output_dir/ASVs2OTUs.log 
wait

### Discard singleton OTUs
if [[ $remove_singletons == "TRUE"  ]]; then
    printf "Discarding singletons ... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --xsize --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error

    mv $output_dir/OTUs.temp.fasta tempdir/
else
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 1 \
    --sizein --xsize --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error

    mv $output_dir/OTUs.temp.fasta tempdir/
fi



#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

#Delete tempdirs if debugger != true
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir ]]; then
        rm -rf tempdir
    fi
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
fi

OTU_count=$(grep -c "^>" $output_dir/OTUs.fasta)
input_ASV_count=$(grep -c "^>" $ASV_fasta_size)
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file
printf "# ASVs clustered to OTUs with vsearch (see 'Core command' below for the used settings).

Clustering formed $OTU_count OTUs 
  [input contained $input_ASV_count ASVs (sequences)].

Files in 'ASVs2OTUs' directory:
# OTUs.fasta            = FASTA formated representative OTU sequences. 
# OTU_table.txt         = OTU distribution table per sample (tab delimited file).
# OTUs.uc               = uclust-like formatted clustering results for OTUs.
# $fastasize.size.fasta = size annotated input sequences [size annotation based on the input table $ASV_table] 

Core command -> 
vsearch $seqsort input.fasta $id $simtype $strands $mask $centroid_in $maxaccepts $cores $otutype OTUs.fasta --fasta_width 0 --sizein
[before clustering, ASV size annotations were fetched from the ASV table]

Total run time was $runtime sec.\n\n
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.23.0 for clustering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
##########################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

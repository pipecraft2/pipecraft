#!/bin/bash

# Chimera filtering with DADA2 removeBimeraDenovo function for ASVs workflow [MIXED orient amplicons].

##########################################################
###Third-party applications:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
##########################################################

#load env variables
readType=${readType}
extension=${fileFormat}
dataFormat=${dataFormat}
workingDir=${workingDir}
echo ${dada2mode}

#load variables
method=${method}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dirs
output_dir1=$"/input/chimeraFiltered_out.dada2"
output_dir2=$"/input/ASVs_out.dada2"

### Check the existance of previous steps, and skip if present



#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Process samples with dada2 removeBimeraDenovo function in R

printf "# Running DADA2 removeBimeraDenovo \n"
Rlog=$(Rscript /scripts/submodules/dada2_chimeraFilt_mixed_wf.R 2>&1)
echo $Rlog > $output_dir1/dada2_chimeraFilt.log 
wait

# Add "ASV" as a 1st col name
if [[ -s $output_dir2/ASVs_table.txt ]]; then
    sed -i "1 s|^|ASV|" $output_dir2/ASVs_table.txt
fi

### NOTE: dada2mode = MIXED -> does not paste out chimeric ASVs per-sample

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    rm $output_dir1/dada2_chimeraFilt.log 
fi
end=$(date +%s)
runtime=$((end-start))

#Make README.txt file (chimeraFiltered_out.dada2)
printf "# Chimeras were filtered out with DADA2 removeBimeraDenovo function (method = $method).

Files in 'chimeraFiltered_out.dada2':
# *.chimFilt_ASVs.fasta = chimera filtered ASVs per sample. 'Size' denotes the abundance of the ASV sequence  
# seq_count_summary.csv = summary of sequence and ASV counts per sample

Total run time was $runtime sec.
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
########################################################" > $output_dir1/README.txt

#Make README.txt file (ASVs_out.dada2)
#count ASVs
ASV_count=$(grep -c "^>" $output_dir2/ASVs.fasta)
nSeqs=$(awk 'BEGIN{FS=OFS="\t"}NR>1{for(i=2;i<=NF;i++) t+=$i; print t; t=0}' $output_dir2/ASVs_table.txt  | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
nCols=$(awk -F'\t' '{print NF; exit}' $output_dir2/ASVs_table.txt)
nSample=$(awk -v NUM=$nCols 'BEGIN {print (NUM-2)}') # -2 cuz 1st column is ASV_ID and 2nd is Sequence

printf "# ASV table was constructed with DADA2 makeSequenceTable function.

Number of ASVs                       = $ASV_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample

Files in 'ASVs_out.dada2' directory:
# ASVs_table.txt                  = ASV-by-sample table (tab delimited file). [denoised and chimera filtered]
# ASVs.fasta                      = FASTA formated representative ASV sequences [denoised and chimera filtered]
# ASVs_table.denoised.nochim.rds  = rds formatted denoised and chimera filtered ASV table (for DADA2)
# ASVs_table.denoised.rds         = rds formatted denoised ASV table (for DADA2)

##################################################################
###Third-party applications for this process [PLEASE CITE]:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
########################################################" > $output_dir2/README.txt

#Done
printf "\nDONE\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir2"
echo "fileFormat=fasta"
echo "readType=single_end"



#!/bin/bash

# Chimera filtering with DADA2 removeBimeraDenovo function for ASVs workflow.

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

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Process samples with dada2 removeBimeraDenovo function in R
printf "# Running DADA2 removeBimeraDenovo \n"
Rlog=$(Rscript /scripts/submodules/dada2_chimeraFilt_wf.R 2>&1)
echo $Rlog > $output_dir1/dada2_chimeraFilt.log 
wait
printf "\n DADA2 removeBimeraDenovo completed \n"

# Add "ASV" as a 1st col name
if [[ -s $output_dir2/ASVs_table.txt ]]; then
    sed -i "1 s|^|ASV|" $output_dir2/ASVs_table.txt
fi


###Compare 'chimera filtered fasta files per sample' and 'NOT chimera filtered fasta files per sample' to paste out only chimeric sequences per sample
echo "pasting chimeric seqs"
#path to denoised-assembled fastas
#path_denoised=$"/input/denoised_assembled.dada2"
path_denoised=$"$workingDir"
printf "\n path_denoised = $path_denoised \n"

#path to chimera filtered fastas
path_chim_filt=$"/input/chimeraFiltered_out.dada2"

#make dir for chimeras.fasta
mkdir $path_chim_filt/chimeras
#make seqs_count_summary.txt
touch $path_chim_filt/chimeras/seq_count_summary.txt
printf "File\tReads\n" > $path_chim_filt/chimeras/seq_count_summary.txt

for chim_filt_file in $path_chim_filt/*.chimFilt_ASVs.fasta; do
    samp_name=$(basename $chim_filt_file | awk 'BEGIN{FS=".chimFilt_ASVs.fasta"} {print $1}')
    corresponding_denoised_file=$(ls $path_denoised | grep "$samp_name")

    #seqkit paste chimeras
    seqkit grep -w 0 -svf \
    <(seqkit seq -s -w 0 $chim_filt_file) \
    <(seqkit seq -w 0 $path_denoised/$corresponding_denoised_file) \
    > $path_chim_filt/chimeras/$samp_name.chimeras.fasta

    #delete if chimeras file is empty
    if [[ ! -s $path_chim_filt/chimeras/$samp_name.chimeras.fasta ]]; then
        printf "$samp_name.chimeras.fasta\t0\n" >> $path_chim_filt/chimeras/seq_count_summary.txt
        rm $path_chim_filt/chimeras/$samp_name.chimeras.fasta
    else
        #count chimeric seqs
        seq_count=$(grep -c "^>" $path_chim_filt/chimeras/$samp_name.chimeras.fasta)
        printf "$samp_name.chimeras.fasta\t$seq_count\n" >> $path_chim_filt/chimeras/seq_count_summary.txt
    fi
done

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
fi

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file (chimeraFiltered_out.dada2)
printf "# Chimera filtering with dada2 removeBimeraDenovo function.

Files in 'chimeraFiltered_out.dada2':
# *.chimFilt_ASVs.fasta   = chimera filtered ASVs per sample. 'Size' denotes the abundance of the ASV sequence  
# ASVs_table.denoised.rds = rds formatted denoised ASV table (for DADA2)
# seq_count_summary.txt   = summary of sequence and ASV counts per sample

Files in 'chimeraFiltered_out.dada2/chimeras':
# *.chimeras.fasta        = ASVs per sample that were flagged as chimeras (and thus discarded).

Core command -> 
removeBimeraDenovo(ASV_tab, method = $method, multithread = FALSE, verbose = FALSE)

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
printf "# Made an ASV table with dada2 makeSequenceTable function.

Number of ASVs = $ASV_count

Files in 'ASVs_out.dada2' directory:
# ASVs_table.txt                  = ASV distribution table per sample (tab delimited file)
# ASVs.fasta                      = FASTA formated representative ASV sequences
# ASVs_table.denoised.nochim.rds  = rds formatted denoised and chimera filtered ASV table (for DADA2)

Core command -> 
makeSequenceTable(merged_paired_end_inputs) [chimeras removed with dada2 removeBimeraDenovo function; see 'chimeraFiltered_out.dada2' directory]

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



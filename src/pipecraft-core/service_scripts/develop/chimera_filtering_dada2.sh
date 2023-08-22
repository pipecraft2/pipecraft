#!/bin/bash

# Chimera filtering with DADA2 removeBimeraDenovo function.
# Input = fastq or fasta files (if fasta, then make dummy fastq files for DADA2 derepFastq)

##########################################################
###Third-party applications:
#dada2 v1.28
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583. https://doi.org/10.1038/nmeth.3869
    #Copyright (C) 2007 Free Software Foundation, Inc.
    #Distributed under the GNU LESSER GENERAL PUBLIC LICENSE
    #https://github.com/benjjneb/dada2
#seqkit v2.3.0
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #Distributed under the MIT License
    #Copyright © 2016-2019 Wei Shen, 2019 Oxford Nanopore Technologies.
    #https://bioinf.shenwei.me/seqkit/
##########################################################


#load env variables
extension=${fileFormat}
workingDir=${workingDir}

#load variables
method=${method}


# If input is fasta, then make dummy fastq for DADA2 'derepFastq'


### Process samples with dada2 removeBimeraDenovo function in R
printf "# Running DADA2 removeBimeraDenovo \n"
errormessage=$(Rscript /scripts/submodules/dada2_chimeraFilt.R 2>&1)
echo $errormessage > $output_dir/R_run.log 
wait
printf "\n DADA2 removeBimeraDenovo completed \n"






#Compare 'chimera filtered fasta files per sample' and 'NOT chimera filtered fasta files per sample' to paste out only chimeric sequences per sample

echo "dada2 bash module for pasting chimeric seqs"

#path to denoised-assembled fastas
path_denoised=$"/input/denoised_assembled.dada2"
#path to chimera filtered fastas
path_chim_filt=$"/input/chimeraFiltered_out.dada2"

echo $path_denoised
echo $path_chim_filt

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



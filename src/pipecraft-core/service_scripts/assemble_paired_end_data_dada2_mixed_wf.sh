#!/bin/bash

# denoise and assemble paired-end data with DADA2 dada and mergePairs functions. For DADA2 full workflow.

##########################################################
###Third-party applications:
#dada2 v1.26
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

#load variables
minOverlap=${minOverlap}
maxMismatch=${maxMismatch}
trimOverhang=${trimOverhang}
justConcatenate=${justConcatenate}
pool=${pool}
qualityType=${qualityType}

#Source for functions
source /scripts/submodules/framework.functions.sh

start=$(date +%s)
for folder in /input/primersCut_out/fwd_orient/qualFiltered_out /input/primersCut_out/rev_orient/qualFiltered_out; do
    ### Check that at least 2 samples are provided
    files=$(ls $folder | grep -c ".$extension")
    if (( $files < 4 )); then
        printf '%s\n' "ERROR]: please provide at least 2 samples for the ASVs workflow
    >Quitting" >&2
        end_process
    fi
    
    #output_dir
    x=$(echo $folder | awk 'BEGIN{FS=OFS="/"}(NF=NF-1) 1')
    output_dir=$"$x/denoised_assembled.dada2"
    export output_dir

    #############################
    ### Start of the workflow ###
    #############################
    workingDir=$folder # workingDir as fwd_orient or rev_orinet
    
    ### Process samples with dada2 removeBimeraDenovo function in R
    printf "# Running DADA2 denoising and assembling \n"
    Rlog=$(Rscript /scripts/submodules/dada2_denoise_assemble_wf.R 2>&1)
    echo $Rlog >> $output_dir/denoise_assemble.log 
    wait
    echo "workingDir=$output_dir"

    #####################################
    ### CLEAN AND COMPILE README FILE ###
    #####################################
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file 
    printf "# Denoising and assembling of PAIRED-END sequencing data with dada2.

    ### NOTE: ### 
    Input sequences must be made up only of A/C/G/T for denoising (i.e maxN must = 0 in quality filtering step). Otherwise DADA2 fails, and no output is generated.
    #############

    Files in 'denoised_assembled.dada2':
    # *.ASVs.fasta   = denoised and assembled ASVs per sample. 'Size' denotes the abundance of the ASV sequence.  
    # Error_rates_R1.pdf    = plots for estimated R1 error rates
    # Error_rates_R2.pdf    = plots for estimated R2 error rates
    # seq_count_summary.csv = summary of sequence and ASV counts per sample
    # *.rds = R objects for dada2.

    Core commands -> 
    learn errors: err = learnErrors(input)
    dereplicate:  derep = derepFastq(input)
    denoise:      dadaFs = dada(input, err = err, pool = pool)
    assemble:     mergePairs(inputR1, dereplicatedR1, inputR2, dereplicatedR2, maxMismatch = $maxMismatch, minOverlap = $minOverlap, justConcatenate = $justConcatenate, trimOverhang = $trimOverhang)

    Total run time was $runtime sec.
    ##################################################################
    ###Third-party applications for this process [PLEASE CITE]:
    #dada2 v1.26
        #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
        #https://github.com/benjjneb/dada2
    ########################################################" > $output_dir/README.txt
done 

#Done
printf "\nDONE\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
if [[ -z $pool ]]; then
    echo "fileFormat=fasta"
else
    echo "fileFormat=$extension"
fi
if [[ -z $pool ]]; then
    echo "readType=single_end"
else
    echo "readType=paired_end"
fi
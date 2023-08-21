#!/bin/bash

#Input = single-end fasta/fastq files.
#Output = FASTA formated zOTU sequences and zOTU_table.txt, and optionally OTU sequences and OTU_table.txt

# Sequence denoising and clustering

##########################################################
###Third-party applications:
#vsearch v2.23.0
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
    #Distributed under the GNU General Public License version 3 by the Free Software Foundation
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422
    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
    #Copyright (C) 2007-2021 Ole Tange, http://ole.tange.dk and Free Software Foundation, Inc.
    #Distributed under the License GPLv3+
#pigz
##########################################################

###############################
###############################
#load variables
id=$"--id ${similarity_threshold}"              # positive float (0-1)  # OTU clustering if id < 1
id_float=${similarity_threshold}
zid=$"--id ${zOTUs_similarity_threshold}"       # positive float (0-1)  # for zOTU table
strands=$"--strand ${strands}"                  # both/plus
minsize=$"--minsize ${minsize}"                 # positive integer (default, 8)

#additional options
unoise_alpha=$"--unoise_alpha ${unoise_alpha}"  # positive integer (default, 2)
denoise_level=${denoise_level}                  # list: "global" or "individual"
chimerarm=${remove_chimeras}                    # TRUE or undefined
cores=$"--threads ${cores}"                     # positive integer
abskew=$"--abskew ${abskew}"                    # positive integer (default, 16)
simtype=$"--iddef ${similarity_type}"           # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
maxaccepts=$"--maxaccepts ${maxaccepts}"        # positive integer (default, 1)
maxrejects=$"--maxrejects ${maxrejects}"        # positive integer (default, 32)
mask=$"--qmask ${mask}"                         # list: --qmask dust, --qmask none
###############################
###############################

#############################
### Start of the workflow ###
#############################
#output dir
output_dir=$"/input/clustering_out"

## Number of cores for GNU parallel. [NOT WORKING]
#NCORES=$cores

start=$(date +%s)
# Source for functions
source /scripts/submodules/framework.functions.sh

### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check single-end data
prepare_SE_env

### Pre-process samples
printf "Checking files ... \n"
for file in *.$fileFormat; do
    #Read file name; without extension
    input=$(echo $file | sed -e "s/.$fileFormat//")
    #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
    check_gz_zip_SE
    ### Check input formats
    check_extension_fastx
done

#If input is FASTQ then convert to FASTA
if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
    for file in *.$extension; do
        samp_name=$(basename $file | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
        checkerror=$(seqkit fq2fa -t dna --line-width 0 $file -o $samp_name.fasta 2>&1)
        check_app_error
    done

    was_fastq=$"true"
    export was_fastq
    extension=$"fasta"
    export extension
fi

#tempdir
if [[ -d tempdir ]]; then
    rm -rf tempdir
fi
mkdir -p tempdir

printf "Dereplication of individual samples ... \n"
## Dereplication of individual samples, add sample ID to the header
derep_rename () {
  samp_name=$(basename $1 | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
  vsearch \
    --derep_fulllength "$1" \
    --relabel_sha1 \
    --output - \
    --fasta_width 0 \
    --sizein --sizeout \
  | sed 's/>.*/&;sample='"$samp_name"'/' > tempdir/"$samp_name".fasta
}
export -f derep_rename
find . -maxdepth 1 -name "*.$extension" | parallel -j 1 "derep_rename {}"

### Global dereplication
printf "Dereplicating globally ... \n"
find tempdir -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
| vsearch \
--derep_fulllength - \
--uc tempdir/Glob_derep.uc \
--output - \
--fasta_width 0 \
--threads 1 \
--sizein --sizeout > $output_dir/Glob_derep.fasta

## Denoizing sequences globally
if [[ $denoise_level == "global" ]]; then
  printf "Denoizing sequences globally ... \n"
 
  ### UNOISE3
  printf "Unoise3 ... \n"
  checkerror=$(vsearch \
  --cluster_unoise $output_dir/Glob_derep.fasta \
  $strands \
  $minsize \
  $unoise_alpha \
  $simtype \
  $mask \
  $maxaccepts \
  $maxrejects \
  $cores \
  --centroids $output_dir/zOTUs.fasta \
  --uc $output_dir/zOTUs.uc \
  --fasta_width 0 \
  --sizein --sizeout 2>&1)
  check_app_error
  
  ## Remove chimera
  printf "Remove chimeras ... \n"

  if [[ $chimerarm == "true" ]]; then
    checkerror=$(vsearch \
    --sortbysize $output_dir/zOTUs.fasta \
    --output - \
    | vsearch --uchime3_denovo - \
    $abskew \
    --nonchimeras $output_dir/zOTUs_noChim.temp.fasta \
    --chimeras $output_dir/UNOISE_Chimeras.fasta 2>&1)
    check_app_error
  
    ## Count number of chimeric sequences
    chimeras=$(grep -c "^>" $output_dir/UNOISE_Chimeras.fasta)

    ## Replace zOTUs with chimera-filtered zOTUs
    rm $output_dir/zOTUs.fasta
    checkerror=$(vsearch --fastx_filter $output_dir/zOTUs_noChim.temp.fasta \
    --fasta_width 0 --fastaout $output_dir/zOTUs.fasta 2>&1)
    check_app_error
    rm $output_dir/zOTUs_noChim.temp.fasta
  fi
fi  # end of global denoising


## Denoizing sequences individually for each sample
if [[ $denoise_level == "individual" ]]; then
  mkdir -p tempdir_denoize
  mkdir -p tempdir_chimera

  ## Function to denoise and remove chimera for each sample individually 
  denoise_and_chim () {
    
    samp_name=$(basename $1)

    ## Denoise sample
    checkerror=$(vsearch \
    --cluster_unoise "$1" \
      $strands \
      $minsize \
      $unoise_alpha \
      $simtype \
      $mask \
      $maxaccepts \
      $maxrejects \
      --threads 1 \
      --centroids tempdir_denoize/"$samp_name" \
      --fasta_width 0 \
      --sizein --sizeout 2>&1)
    check_app_error
    
    ## Remove chimera
    if [[ $chimerarm == "true" ]]; then
      checkerror=$(vsearch \
        --sortbysize tempdir_denoize/"$samp_name" \
        --output - \
        | vsearch \
        --uchime3_denovo - \
        $abskew \
        --nonchimeras tempdir_chimera/NonChim_"$samp_name" \
        --chimeras tempdir_chimera/Chim_"$samp_name" 2>&1)
      check_app_error
    fi
  }

  export -f denoise_and_chim
  export -f check_app_error
  export -f end_process

  export chimerarm="$chimerarm"
  export strands="$strands"
  export minsize="$minsize"
  export unoise_alpha="$unoise_alpha"
  export simtype="$simtype"
  export mask="$mask"
  export maxaccepts="$maxaccepts"
  export maxrejects="$maxrejects"

  ## Take dereplicated samples and apply denoising function
  printf "Denoizing sequences individually ... \n"
  find tempdir -maxdepth 1 -name "*.fasta" | parallel -j 1 "denoise_and_chim {}"

  if [[ $chimerarm == "true" ]]; then
    printf "Removing chimeras ... \n"
    find tempdir_chimera -maxdepth 1 -name "Chim_*.fasta" | parallel -j 1 "cat {} >> tempdir/All_chimera.fasta"
    ## Count chimeric sequences
    chimeras=$(grep -c "^>" tempdir/All_chimera.fasta)

    ## Combine and dereplicate denoised sequences
    find tempdir_chimera -maxdepth 1 -name "NonChim_*.fasta" | parallel -j 1 "cat {}" \
    | vsearch \
    --derep_fulllength - \
    --output $output_dir/zOTUs.fasta \
    --uc $output_dir/zOTUs.uc \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout
  
  else
    ## Combine and dereplicate denoised sequences (without chimera removal step)
    find tempdir_denoize -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
    | vsearch \
    --derep_fulllength - \
    --output $output_dir/zOTUs.fasta \
    --uc $output_dir/zOTUs.uc \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout
  fi
  
fi # end of individual denoising

### OTU tables
### Cat dereplicated individual samples for making an OTU table
cat tempdir/*.fasta > $output_dir/Dereplicated_samples.fasta

## Prepare table with sequence abundance per sample
seqkit seq --name $output_dir/Dereplicated_samples.fasta \
  | awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
  | sed 's/size=//; s/sample=//' \
  > tempdir/ASV_table_long.txt

## zOTU table creation
printf "Making zOTU table ... \n"
Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
  --derepuc      tempdir/Glob_derep.uc \
  --uc           "$output_dir"/zOTUs.uc \
  --asv          tempdir/ASV_table_long.txt \
  --rmsingletons FALSE \
  --output       "$output_dir"/zOTU_table.txt 2>&1)
echo $Rlog > tempdir/zOTU_table_creation.log 
wait

## Perform OTU clustering (if required, id < 1)
if [[ $id_float != 1 ]]; then
  printf "\n Clustering zOTUs ... \n"

  ### Clustering
  checkerror=$(vsearch --cluster_size \
  $output_dir/zOTUs.fasta \
  $id \
  $simtype \
  $strands \
  $mask \
  $maxaccepts \
  $maxrejects \
  $cores \
  --centroids $output_dir/OTUs.fasta \
  --uc $output_dir/OTUs.uc \
  --fasta_width 0 \
  --sizein --sizeout 2>&1)
  check_app_error

  ## OTU table creation
  printf "Making OTU table ... \n"
  Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
    --derepuc      tempdir/Glob_derep.uc \
    --uc           "$output_dir"/OTUs.uc \
    --asv          tempdir/ASV_table_long.txt \
    --rmsingletons FALSE \
    --output       "$output_dir"/OTU_table.txt 2>&1)
  echo $Rlog > tempdir/OTU_table_creation.log 
  wait
fi # end of OTU clustering

### remove ";sample=.*;" from OTU.fasta files
if [[ -f $output_dir/OTUs.fasta ]]; then
  sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
fi
if [[ -f $output_dir/zOTUs.fasta ]]; then
  sed -i 's/;sample=.*;/;/' $output_dir/zOTUs.fasta
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

mv $output_dir/Glob_derep.fasta tempdir/Glob_derep.fasta
mv $output_dir/Dereplicated_samples.fasta tempdir/Dereplicated_samples.fasta

#If input = FASTQ, then mkdir for converted fasta files
if [[ $was_fastq == "true" ]]; then
    mkdir -p $output_dir/clustering_input_to_FASTA
    mv *.fasta $output_dir/clustering_input_to_FASTA
fi

#Delete decompressed files if original set of files were compressed
if [[ $check_compress == "gz" ]] || [[ $check_compress == "zip" ]]; then
  delete_uncompressed=$(echo $fileFormat | (awk 'BEGIN{FS=OFS="."} {print $(NF-1)}';))
  rm *.$delete_uncompressed
fi

#Delete tempdirs
if [[ $debugger != "true" ]]; then
  if [[ -d tempdir ]]; then
      rm -rf tempdir
  fi
  if [[ -d tempdir2 ]]; then
      rm -rf tempdir2
  fi
  if [[ -d tempdir_denoize ]]; then
      rm -rf tempdir_denoize
  fi
  if [[ -d tempdir_chimera ]]; then
      rm -rf tempdir_chimera
  fi
  if [[ -f $output_dir/zOTU_table_creation.log ]]; then
      rm -f $output_dir/zOTU_table_creation.log
  fi
  if [[ -f $output_dir/OTU_table_creation.log ]]; then
      rm -f $output_dir/OTU_table_creation.log
  fi
else 
  #compress files in /tempdir
  pigz tempdir/*
fi

#Make README.txt file
size_zotu=$(grep -c "^>" $output_dir/zOTUs.fasta)
end=$(date +%s)
runtime=$((end-start))

printf "Sequence denoising formed $size_zotu zOTUs (zero-radius OTUs).

Files in 'clustering_out' directory:
# zOTUs.fasta    = FASTA formated denoized sequences (zOTUs.fasta). Headers are renamed according to sha1 algorithm in vsearch.
# zOTU_table.txt = zOTU distribution table per sample (per input file in the working directory).
# zOTUs.uc       = uclust-like formatted clustering results for zOTUs

Core command -> 
UNOISE: vsearch --cluster_unoise dereplicated_sequences.fasta $strands $minsize $unoise_alpha $simtype $mask $maxaccepts $maxrejects $cores --centroids zOTUs.fasta --uc zOTUs.uc \n\n" > $output_dir/README.txt

## If additional clustering was performed
if [[ $id_float != 1 ]]; then
    size_otu=$(grep -c "^>" $output_dir/OTUs.fasta)
    printf "Additional clustering of zOTUs at $id similarity threshold formed $size_otu OTUs.
    # OTUs.fasta    = FASTA formated representative OTU sequences. Headers are renamed according to sha1 algorithm in vsearch.
    # OTU_table.txt = OTU distribution table per sample (per input file in the working directory).
    # OTUs.uc       = uclust-like formatted clustering results for OTUs.
    
    Core command -> 
    clustering: vsearch --cluster_size zOTUs.fasta $id $simtype $strands $mask $maxaccepts $maxrejects $cores --centroids OTUs.fasta --uc OTUs.uc \n\n" >> $output_dir/README.txt
fi

## Chimera stats
if [[ $chimerarm == "true" ]]; then
    printf "Chimera removal step eliminated $chimeras sequences\n" >> $output_dir/README.txt
fi

## if input was fastq
if [[ $was_fastq == "true" ]]; then
  printf "\nInput was fastq; converted those to fasta before clustering. 
  Converted fasta files in directory 'clustering_input_to_FASTA' \n" >> $output_dir/README.txt
fi

printf "\nIf samples are denoised individually rather by pooling all samples together, 
reducing minsize to 4 is more reasonable for higher sensitivity.
\n" >> $output_dir/README.txt

printf "\nTotal run time was $runtime sec.\n\n
##################################################################
###Third-party applications for this process [PLEASE CITE]:
#vsearch v2.23.0 for clustering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422 for job parallelisation 
    #Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
##########################################################" >> $output_dir/README.txt

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single-end"

#!/bin/bash

## DEICODE (Robust Aitchison PCA on sparse compositional metabarcoding data)
# Step 1. rCLR
#   Data transformation using centered log ratio on only non-zero values (no pseudo count added)
# Step 2. RPCA
#   Dimensionality reduction through robust PCA on only the non-zero values of the data
#   Sparse data handled through the use of matrix completion (all zeros are treated as missing values)

## See method description in Martino et al., 2019 (DOI 10.1128/mSystems.00016-19)

##########################################################
###Third-party applications:
#DEICODE v0.2.4
    #citation: Martino C., Morton J.T., Marotz C.A., Thompson L.R., Tripathi A., Knight R., et al. (2019). A Novel Sparse Compositional Technique Reveals Microbial Perturbations. mSystems. doi: 10.1128/mSystems.00016-19. 
    #Distributed under the BSD-3 license
    #https://github.com/biocore/DEICODE
#biom-format v2.1.12
    #citation: McDonald D., Clemente J.C., Kuczynski J., Rideout J.R., Stombaugh J., Wendel D., et al. (2012). The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 1, 2047-217X-1–7. doi: 10.1186/2047-217X-1-7.
    #Distributed under the Modified BSD License
    #https://biom-format.org/
##################################################################

## Input:
# otu_table = tab-delimited OTU table (columns = samples, rows = OTUs)
# otu_subset = optional, (one ID per row)
# min_otu_reads = threshold for the minimum number of reads per OTU across all samples
# min_sample_reads = threshold for the minimum number of reads per sample across all OTUs

## Dependencies:
# subset_otu_table.R
# rpca_script.py

#env variables
workingDir=${workingDir}
extension=$fileFormat
## Load variables
subset_IDs=${subset_IDs}
min_otu_reads=${min_otu_reads}
min_sample_reads=${min_sample_reads}

#Source for functions
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/DEICODE_out"

#Automatic search for OTU_table.txt or ASVs_table.txt (standard PipeCraft2 output file names), otherwise use the file that was specified in the panel
if [[ -e "$workingDir/OTU_table.txt" ]]; then
    otu_table=$"$workingDir/OTU_table.txt"
    printf "\n input table = $otu_table \n"
elif [[ -e "$workingDir/ASVs_table.txt" ]]; then
    otu_table=$"$workingDir/ASVs_table.txt"
    printf "\n input table = $otu_table \n"
elif [[ $table == "undefined" ]]; then
    printf '%s\n' "ERROR]: input table was not specified and cannot find OTU_table.txt or ASVs_table.txt in the working dir.
    >Quitting" >&2
    end_process
else
    #get input OTU table file
    regex='[^/]*$'
    otu_table_temp=$(echo $table | grep -oP "$regex")
    otu_table=$(printf "/extraFiles/$otu_table_temp")
    printf "\n input table = $otu_table \n "
fi

#If specified, get OTUs subset
if [[ $subset_IDs == "undefined" ]]; then
    :
else
    regex='[^/]*$'
    subset_IDs_temp=$(echo $subset_IDs | grep -oP "$regex")
    otu_subset=$(printf "/extraFiles2/$subset_IDs_temp")
    printf "\n subset = $otu_subset \n "
fi

## Outputs
output_biom="$output_dir"/otutab.biom
output_deicode_full="$output_dir"/full
output_deicode_subs="$output_dir"/subs
otu_table_subs="$output_dir"/rclr_subset.tsv
#Remove 'old' output_dir if exist and make new empty one
if [ -d $output_dir ]; then
    rm -rf $output_dir
fi
mkdir -p "$output_dir"

#start time
start=$(date +%s)

## Prepare BIOM file
printf "# Preparing BIOM file\n "
checkerror=$(biom convert \
            -i "$otu_table" \
            -o "$output_biom" \
            --table-type="OTU table" \
            --to-hdf5 2>&1)
check_app_error

## rCLR and RPCA on an OTU table with all OTUs
printf "# Running DEICODE on a full OTU table\n "
checkerror=$(deicode auto-rpca \
                    --in-biom "$output_biom" \
                    --min-feature-count "$min_otu_reads" \
                    --min-sample-count "$min_sample_reads" \
                    --output-dir "$output_deicode_full" 2>&1)
check_app_error

# rCLR-transformed table
output_rclr="$output_deicode_full"/rclr.tsv


## RPCA on a subset of OTU table
if [[ $subset_IDs != "undefined" ]]; then
    printf "# Preparing a subset of rCLR-transformed OTU table\n "

    Rscript /scripts/submodules/subset_otu_table.R \
      --input "$output_rclr" \
      --subset "$otu_subset" \
      --output "$otu_table_subs"

    printf "# Running DEICODE on a subset of OTU table\n "
    python /scripts/submodules/rpca_script.py "$otu_table_subs" "$output_deicode_subs"

else 
    printf "# No OTU IDs specified\n "
fi

#DONE
end=$(date +%s)
runtime=$((end-start))

###Make README.txt file
printf "## DEICODE (Robust Aitchison PCA on sparse compositional metabarcoding data)

# Step 1. rCLR
#   Data transformation using centered log ratio on only non-zero values (no pseudo count added)
# Step 2. RPCA
#   Dimensionality reduction through robust PCA on only the non-zero values of the data
#   Sparse data handled through the use of matrix completion (all zeros are treated as missing values)

## See method description in Martino et al., 2019 (DOI 10.1128/mSystems.00016-19)

Files in 'DEICODE_out' directory:
#   - otutab.biom          =  full OTU table in BIOM format
#   - rclr_subset.tsv      =  rCLR-transformed subset of OTU table *
# DEICODE_out/full/
#   - distance-matrix.tsv  =  distance matrix between the samples, based on full OTU table
#   - ordination.txt       =  ordination scores for samples and OTUs, based on full OTU table
#   - rclr.tsv             =  rCLR-transformed OTU table
# DEICODE_out/subs/
#   - distance-matrix.tsv  =  distance matrix between the samples, based on a subset of OTU table *
#   - ordination.txt       =  ordination scores for samples and OTUs, based a subset of OTU table *
#
# *, files are present only if 'subset_IDs' variable was specified

Core commands -> 
biom convert -i $otu_table -o $output_biom --table-type='OTU table' --to-hdf5
deicode auto-rpca --in-biom $output_biom --min-feature-count $min_otu_reads --min-sample-count $min_sample_reads --output-dir $output_deicode_full

Total run time was $runtime sec.

##########################################################
###Third-party applications [PLEASE CITE]:
#DEICODE v0.2.4
    #citation: Martino C., Morton J.T., Marotz C.A., Thompson L.R., Tripathi A., Knight R., et al. (2019). A Novel Sparse Compositional Technique Reveals Microbial Perturbations. mSystems. doi: 10.1128/mSystems.00016-19. 
#biom-format v2.1.12 (for preparing BIOM file)
    #citation: McDonald D., Clemente J.C., Kuczynski J., Rideout J.R., Stombaugh J., Wendel D., et al. (2012). The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 1, 2047-217X-1–7. doi: 10.1186/2047-217X-1-7.
##################################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
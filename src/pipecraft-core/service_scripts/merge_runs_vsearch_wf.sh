#!/bin/bash

# Merge sequencing runs processed with vsearch OTUs workflow if working with multuple runs in multiRunDir. 
 # vsearch clustering for all fasta files in either chimeraFiltered_out or ITSx_out directories.
  # + apply 'curate otu table' when this is applied.
 # Samples with the same name across runs are not automatically merged together; each sample will be tagged with RunID__SampleID.

 # 1. clustering all samples from all runs.
 # 2. Split tables per run
 # 3. curate tables (tj + lenFilt) if 'curate otu table' is enabled
 # 4. merge OTU tables

################################################
###Third-party applications:
# vsearch
################################################
# Checking tool versions
printf "# Checking tool versions ...\n"
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# vsearch version: $vsearch_version\n"

start_time=$(date)
start=$(date +%s)
# source for functions
source /scripts/submodules/framework.functions.sh

# Excecute only if multiDir = true
if [[ ! -d "/input/multiRunDir" ]]; then
    printf '%s\n' "ERROR]: multiRunDir not detected. Cannot merge sequencing runs. 
    >DONE." >&2
    end_process
elif [[ $merge_runs == "true" ]]; then
    printf "Starting merge sequencing runs...\n"
    #output dir
    output_dir=$"/input/multiRunDir/merged_runs"
    export output_dir
    # remove output dir if it already exists
    if [[ -d "$output_dir" ]]; then
        rm -rf $output_dir
    fi
    # create new output dir
    mkdir -p $output_dir
    echo "input tables: $output_feature_table" 
    echo "input fasta: $output_fasta"
    echo "output dir: $output_dir"
else
    printf '%s\n' "ERROR]: Merge sequencing runs is not enabled. Exiting.\n" >&2
    end_process
fi


### 1. clustering all samples from all runs.
 # using dereplicated individual samples from clustering step (chimeraFiltered_out/dereplicated_sequences/*.fasta)

# get run names
if [[ -d "/input/multiRunDir" ]]; then
    cd /input/multiRunDir
    RUNS=$(find . -maxdepth 1 -mindepth 1 -type d | grep -v "tempdir" | grep -v "skip_" | grep -v "merged_runs" | sed -e "s/^\.\///")
    echo "Sequencing runs:"
    echo $RUNS
fi

# add run name to the sample names
printf "Adding run name to the sample names ... \n"
mkdir -p /input/multiRunDir/merged_runs/dereplicated_sequences
for seqrun in $RUNS; do
    mkdir -p /input/multiRunDir/merged_runs/${seqrun}_dereplicated_sequences
    cp /input/multiRunDir/$seqrun/chimeraFiltered_out/dereplicated_sequences/*.fasta /input/multiRunDir/merged_runs/${seqrun}_dereplicated_sequences/
    for file in /input/multiRunDir/merged_runs/${seqrun}_dereplicated_sequences/*.fasta; do
        sed -i "s/sample=/sample=${seqrun}__/g" "$file"
    done
    mv /input/multiRunDir/merged_runs/${seqrun}_dereplicated_sequences/ /input/multiRunDir/merged_runs/
    rm -rf /input/multiRunDir/merged_runs/${seqrun}_dereplicated_sequences
done

### Global dereplication
if [[ -d $output_dir/tempdir ]]; then
    rm -rf $output_dir/tempdir
fi
mkdir -p $output_dir/tempdir

printf "Dereplicating globally ... \n"
find $output_dir/dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
| vsearch \
    --derep_fulllength - \
    --output $output_dir/Glob_derep.fasta \
    --uc $output_dir/tempdir/Glob_derep.uc \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout

### Clustering -- POOLELI - siin needs to pass all CLUSTERING properties. Write to a file? !
printf "Clustering ... \n"
checkerror=$(vsearch $seqsort \
    $output_dir/Glob_derep.fasta \
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

### Cat dereplicated individual samples for making an OTU table
cat dereplicated_sequences/*.fasta > tempdir/Dereplicated_samples.fasta

## Prepare table with sequence abundance per sample
seqkit seq --name tempdir/Dereplicated_samples.fasta \
| awk -F ";" '{print $3 "\t" $1 "\t" $2}' \
| sed 's/size=//; s/sample=//' \
> tempdir/ASV_table_long.txt

### OTU table creation
printf "Making OTU table ... \n"
Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
--derepuc      tempdir/Glob_derep.uc \
--uc           "$output_dir"/OTUs.uc \
--asv          tempdir/ASV_table_long.txt \
--rmsingletons $remove_singletons \
--output       "$output_dir"/OTU_table.txt 2>&1)
echo $Rlog > tempdir/OTU_table_creation.log 
wait

### Discard singleton OTUs
if [[ $remove_singletons == "true"  ]]; then
    printf "Discarding singletons ... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --sizeout --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error
    # remove ";sample=.*;" from OTU.fasta file.
    sed -i 's/;sample=.*;/;/' $output_dir/OTUs.fasta
    # removing ";size=" because OTU table does not have "size" annotations; so the files would fit to LULU
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta 
    mv $output_dir/OTUs.temp.fasta tempdir/
else
    sed -e 's/;sample=.*;/;/' $output_dir/OTUs.temp.fasta > $output_dir/OTUs.fasta
    sed -i 's/;size=.*//' $output_dir/OTUs.fasta
    mv $output_dir/OTUs.temp.fasta tempdir/
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

mv $output_dir/Glob_derep.fasta tempdir/Glob_derep.fasta

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
if [[ $multiDir != "TRUE" ]]; then
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir ]]; then
            rm -rf tempdir
        fi
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        if [[ -f $output_dir/OTU_table_creation.log ]]; then
            rm -f $output_dir/OTU_table_creation.log
        fi
        else 
        compress files in /tempdir
        if [[ -d tempdir ]]; then
            pigz tempdir/*
        fi
    fi
fi

#Make README.txt file
count_features "$output_dir/OTU_table.txt"

end=$(date +%s)
runtime=$((end-start))

printf "# Reads were clustered to OTUs using vsearch (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'clustering_out' directory:
-----------------------------------
# OTUs.fasta    = FASTA formated representative OTU sequences. OTU headers are renamed according to sha1 algorithm in vsearch.
# OTU_table.txt = OTU distribution table per sample (tab delimited file). OTU headers are renamed according to sha1 algorithm in vsearch.
# OTUs.uc       = uclust-like formatted clustering results for OTUs.

Number of OTUs                       = $feature_count
Number of sequences in the OTU table = $nSeqs
Number of samples in the OTU table   = $nSample

Core command -> 
clustering: vsearch $seqsort dereplicated_sequences.fasta $id $simtype $strands $mask $centroid_in $maxaccepts $cores $otutype OTUs.fasta " > $output_dir/README.txt

## if input was fastq
if [[ $was_fastq == "true" ]]; then
printf "\n\nInput was fastq; converted those to fasta before clustering. 
Converted fasta files in directory 'clustering_input_to_FASTA' \n" >> $output_dir/README.txt
fi

printf "

##############################################
###Third-party applications for this process:
#vsearch (version $vsearch_version)
#citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
#https://github.com/torognes/vsearch
#GNU Parallel 20210422 for job parallelisation 
#Citation: Tange, O. (2021, April 22). GNU Parallel 20210422 ('Ever Given'). Zenodo. https://doi.org/10.5281/zenodo.4710607
################################################" >> $output_dir/README.txt

















# check if 2nd column of the OTU table is "Sequence"
if [[ $(awk -F'\t' 'NR==1 {print ($2=="Sequence")}' "$output_feature_table") == 1 ]]; then
    echo "2nd column is sequence, ok"
else
    echo "2nd column is not sequence, adding sequences to the table"
    Rscript /scripts/submodules/add_sequences_to_table.R $output_feature_table $output_fasta
     # output from latter is ${output_feature_table%.*}_wSeqs.txt
    # add sequences to the table
    fastasize=$(basename $output_fasta | awk 'BEGIN{FS=OFS="."}NF{NF -=1}1')
    awk 'NR>1{for(i=3;i<=NF;i++) t+=$i; print ">"$1";size="t"\n"$2; t=0}' ${output_feature_table%.*}_wSeqs.txt > $output_dir/$fastasize.size.fasta
fi

    #export for R
    ASV_fasta_size=$"$output_dir/$fastasize.size.fasta"
    export ASV_fasta_size




















    ### Convert wide format OTU table to long format in awk (bash) 
    #  handles "Sequence" column


    checkerror=$(cat $table_file \
    | awk '
    BEGIN { 
    FS="\t"; OFS="\t";
    print "OTU", "SampleID", "Abundance"   # Header of the resulting table
    }
    NR==1 {
    seq_col = -1;                      # Initialize sequence column index
    for (i=2; i<=NF; i++) {
        if ($i == "Sequence") {        # Check if column is "Sequence"
        seq_col = i;                   # Store sequence column index
        continue;                      # Skip this column
        }
        sampleIDs[i] = $i;             # Store sample IDs from the header row
    }
    }
    NR>1 {
    otu = $1;                          # Get the OTU ID from the first column
    for (i=2; i<=NF; i++) {
        if (i == seq_col) continue;    # Skip the sequence column
        if ($i > 0) {                  # Skip zero abundances
        print otu, sampleIDs[i], $i;   # Print OTU ID, SampleID, and Abundance
        }
    }
    }' > $output_dir/table_long.txt 2>&1)
    check_app_error



### Make README.txt file (merged_runs)
end=$(date +%s)
runtime=$((end-start))
printf "# Merged sequencing runs with DADA2 mergeSequenceTables function.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input tables:\n" > $output_dir/README.txt

# Add each input table path
IFS=',' read -ra TABLES <<< "$output_feature_table"
for table in "${TABLES[@]}"; do
    printf "%s\n" "$table" >> $output_dir/README.txt
done

printf "
Output files:
------------
# ASVs_table.txt = merged ASV abundance table
# ASVs.fasta     = merged ASV sequences

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample " >> $output_dir/README.txt


if [[ $collapseNoMismatch == "true" ]]; then
    printf "\n
Outputs after CollapsedNoMismatch = true:
-----------------------------------------
$ASVs_collapsed_result \n" >> $output_dir/README.txt
fi

if [[ $collapseNoMismatch == "true" ]] && [[ -f $output_dir/${feature_table_base_name%%.txt}_collapsed.txt ]]; then
    count_features $output_dir/${feature_table_base_name%%.txt}_collapsed.txt
    printf "
${feature_table_base_name%%.txt}_collapsed.txt = merged and collapsed ASV abundance table
${fasta_base_name%%.fasta}_collapsed.fasta = merged and collapsed ASV sequences

Number of ASVs                       = $feature_count
Number of sequences in the ASV table = $nSeqs
Number of samples in the ASV table   = $nSample " >> $output_dir/README.txt
fi

printf "\n
#############################################
###Third-party applications for this process:
# dada2 (version $dada2_version)
    #citation: Callahan, B., McMurdie, P., Rosen, M. et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581-583. https://doi.org/10.1038/nmeth.3869
    #https://github.com/benjjneb/dada2
# vsearch (version $vsearch_version)
    #citation: Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584
    #https://github.com/torognes/vsearch
#############################################" >> $output_dir/README.txt

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

# variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"

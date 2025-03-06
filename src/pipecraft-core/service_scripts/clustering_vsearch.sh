#!/bin/bash

# Sequence clustering with vsearch
#Input = single-end fasta/fastq files.
#Output = FASTA formated representative OTU sequences and OTU_table.txt.

################################################
###Third-party applications:
#vsearch v2.23.0
    #https://github.com/torognes/vsearch
#GNU Parallel 20210422
# seqkit v2.3.0
#pigz v2.4
################################################
# checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"
printf "# pipeline = $pipeline\n"
printf "# service = $service\n"
printf "# clustering dirs = $workingDir\n"

#load variables
id=$"--id ${similarity_threshold}"          # positive float (0-1)
otutype=$"--${OTU_type}"                    # list: --centroids, --consout
strands=$"--strand ${strands}"              # list: --strand both, --strand plus
remove_singletons=$"${remove_singletons}"   # true/false

#additional options
seqsort=$"${sequence_sorting}"           # list: --cluster_size or --cluster_fast, --cluster_smallmem
simtype=$"--iddef ${similarity_type}"    # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
centroid=$centroid_type                  # list: similarity, abundance
maxaccepts=$"--maxaccepts ${maxaccepts}" # pos integer
mask=$"--qmask ${mask}"                  # list: --qmask dust, --qmask none
cores=$"--threads ${cores}"              # pos integer
###############################
# Source for functions
source /scripts/submodules/framework.functions.sh

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

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "vsearch OTUs with multiple sequencing runs in multiRunDir"
    echo "Process = vsearch clustering"
    cd /input/multiRunDir

    # rm old clustering params if they exist
    if [[ -f "/input/multiRunDir/.clustering_params" ]]; then
        rm /input/multiRunDir/.clustering_params
    fi
    if [[ -f "/input/multiRunDir/.derep_seqs_dirs" ]]; then
        rm /input/multiRunDir/.derep_seqs_dirs
    fi

    # check if prev_step.temp exists
    if [[ -f "$workingDir/.prev_step.temp" ]]; then
        prev_step=$(cat $workingDir/.prev_step.temp) # for checking previous step (output temp file from ITS_extractor.sh)
        printf "# prev_step = $prev_step\n"
        ITS_region_for_clustering=$(cat $workingDir/.ITS_region_for_clustering.temp)
        ITS_full_and_partial=$(cat $workingDir/.ITS_full_and_partial.temp)
        printf "# ITS_region_for_clustering = $ITS_region_for_clustering\n"
        printf "# ITS_full_and_partial = $ITS_full_and_partial (=FALSE if empty)\n"
        rm $workingDir/.prev_step.temp
        rm $workingDir/.ITS_region_for_clustering.temp
        rm $workingDir/.ITS_full_and_partial.temp
    fi
    
    # if working with multiRunDir, and the previous step was ITSx and full_and_partial = FALSE
    if [[ "$prev_step" = "ITSx" ]] && [[ "$ITS_full_and_partial" == "" ]]; then 
        DIRS=$(find . -maxdepth 3 -mindepth 1 -type d | grep "ITSx_out/${ITS_region_for_clustering}" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    # if working with multiRunDir, and the previous step was ITSx and full_and_partial = TRUE
    elif [[ "$prev_step" = "ITSx" ]] && [[ "$ITS_full_and_partial" == "full_and_partial" ]]; then   
        DIRS=$(find . -maxdepth 4 -mindepth 1 -type d | grep "ITSx_out/${ITS_region_for_clustering}/${ITS_full_and_partial}" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///")
    # if working with multiRunDir, and the previous step was not ITSx
    else
        DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "chimeraFiltered_out" | grep -v "chimeraFiltered_out.dada2"| grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    fi
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = vsearch clustering"
    DIRS=$(pwd)
    printf "\n workingDir = $DIRS \n"
fi

#############################
### Start of the workflow ###
#############################
# Store output files in arrays for multiple runs
declare -a output_feature_tables
declare -a output_fastas

### looping through multiple sequencing runs (dirs in multiRunDir) if the $WD=multiRunDir, otherwise just doing single seqrun analyses
for seqrun in $DIRS; do
    start_time=$(date)
    start=$(date +%s)

    cd $seqrun
    if [[ $multiDir == "TRUE" ]]; then
        ### Check if the dir has the specified file extension; if not then ERROR
        count=$(ls -1 *.$fileFormat 2>/dev/null | wc -l)
        if [[ $count == 0 ]]; then
            printf '%s\n' "ERROR]: cannot find files with specified extension '$fileFormat' in dir $seqrun.
            Please check the extension of your files and specify again.
            >Quitting" >&2
            end_process
        fi
        #output dir
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/clustering_out"
        export output_dir
        # Store output files in arrays for multiple runs
        output_feature_tables+=("$output_dir/OTU_table.txt")
        output_fastas+=("$output_dir/OTUs.fasta")

        ### Prepare working env and check single-end data
        first_file_check
        prepare_SE_env
    else
        #output dir
        output_dir=$"/input/clustering_out"
        export output_dir
        # Store single output files
        output_feature_table="$output_dir/OTU_table.txt"
        output_fasta="$output_dir/OTUs.fasta"
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check single-end data
        prepare_SE_env
    fi

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

    if [[ ! -d dereplicated_sequences ]]; then
        rm -rf dereplicated_sequences
    fi
    mkdir -p dereplicated_sequences

    ### Rename sequences to sha1
        # and dereplication of individual samples, add sample ID to the header
    derep_rename () {
    samp_name=$(basename $1 | awk 'BEGIN{FS="."} {$NF=""; print $0}' | sed 's/ //g')
    vsearch \
        --derep_fulllength "$1" \
        --relabel_sha1 \
        --output - \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout \
    | sed 's/>.*/&;sample='"$samp_name"'/' > dereplicated_sequences/"$samp_name".fasta
    }
    export -f derep_rename
    printf "Dereplication of individual samples ... \n"
    find . -maxdepth 1 -name "*.$extension" | parallel -j 1 "derep_rename {}"

    ### Global dereplication
    printf "Dereplicating globally ... \n"
    find dereplicated_sequences -maxdepth 1 -name "*.fasta" | parallel -j 1 "cat {}" \
    | vsearch \
    --derep_fulllength - \
    --output $output_dir/Glob_derep.fasta \
    --uc tempdir/Glob_derep.uc \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout

    ### Clustering
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

    ### OTU table creation (--fasta is added to write OTU seqs to OTU_table.txt)
    printf "Making OTU table ... \n"
    Rlog=$(Rscript /scripts/submodules/ASV_OTU_merging_script.R \
    --derepuc      tempdir/Glob_derep.uc \
    --uc           "$output_dir"/OTUs.uc \
    --asv          tempdir/ASV_table_long.txt \
    --rmsingletons $remove_singletons \
    --fasta        $output_dir/OTUs.temp.fasta \
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
    if [[ $debugger != "true" ]]; then
        if [[ -d tempdir ]]; then
            rm -rf tempdir
        fi
        if [[ -d tempdir2 ]]; then
            rm -rf tempdir2
        fi
        else 
        # compress files in /tempdir
        if [[ -d tempdir ]]; then
            pigz tempdir/*
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
    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done

# Validate and combine output files
if [[ $multiDir == "TRUE" ]]; then
    # Create an array of dereplicated_sequences directories (for merge_runs_vsearch_wf.sh)
    derep_dirs=()
    for dir in $DIRS; do
        derep_dir="${dir}/dereplicated_sequences"
        derep_dirs+=("$derep_dir")
    done
    printf "%s\n" "${derep_dirs[@]}" > /input/multiRunDir/.derep_seqs_dirs

    # Check each file in the arrays
    valid_feature_tables=()
    valid_fastas=()
    
    for table in "${output_feature_tables[@]}"; do
        if [[ -f "$table" && -s "$table" ]]; then
            valid_feature_tables+=("$table")
        else
            printf "Warning: Feature table not found or empty: $table\n" >&2
        fi
    done
    
    for fasta in "${output_fastas[@]}"; do
        if [[ -f "$fasta" && -s "$fasta" ]]; then
            valid_fastas+=("$fasta")
        else
            printf "Warning: FASTA file not found or empty: $fasta\n" >&2
        fi
    done
      
    output_feature_table=$(IFS=,; echo "${valid_feature_tables[*]}")
    output_fasta=$(IFS=,; echo "${valid_fastas[*]}")
else
    # Check single run output files
    if [[ ! -f "$output_feature_table" || ! -s "$output_feature_table" ]]; then
        printf "Error: Feature table not found or empty: $output_feature_table\n" >&2
        exit 1
    fi
    if [[ ! -f "$output_fasta" || ! -s "$output_fasta" ]]; then
        printf "Error: FASTA file not found or empty: $output_fasta\n" >&2
        exit 1
    fi
fi

if [[ $multiDir == "TRUE" ]]; then
# write clustering parameters into file for the merge_runs_vsearch_wf.sh
    cat > /input/multiRunDir/.clustering_params << EOF
id="${id}"
otutype="${otutype}"
strands="${strands}"
remove_singletons="${remove_singletons}"
seqsort="${seqsort}"
simtype="${simtype}"
centroid_in="${centroid_in}"
maxaccepts="${maxaccepts}"
mask="${mask}"
cores="${cores}"
EOF
fi

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
echo "output_feature_table=$output_feature_table"
echo "output_fasta=$output_fasta"


#!/bin/bash

# Chimera filtering
# Input = single-end fasta/fastq files. FASTQ files will be converted to FASTA files; output is only FASTA.

################################################
###Third-party applications:
#vsearch v2.23.0
#seqkit v2.3.0
#pigz v2.4
################################################
# checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
seqkit_version=$(seqkit version 2>&1 | awk '{print $2}')
printf "# vsearch (version $vsearch_version)\n"
printf "# seqkit (version $seqkit_version)\n"

#load variables
id=$"--id ${pre_cluster}"                            #float (0-1)
minuniquesize=$"--minuniquesize ${min_unique_size}"  #pos int >0
denovo=${denovo}                                     #false or true
cores=$"--threads ${cores}"                          #pos int
abskew=$"--abskew ${abundance_skew}"                 #pos int
minh=$"--minh ${min_h}"                              #float (0-1)

#Source for functions
source /scripts/submodules/framework.functions.sh

#load path to reference database (if specified)
if [[ $reference_based == "undefined" ]]; then
    :
else
    regex='[^/]*$'
    ref=$(echo $reference_based | grep -oP "$regex")
    db=$(printf "/extraFiles4/$ref")
    database=$db
fi

#ERROR if both chimera filtering methods = false
if [[ $reference_based == "undefined" ]] && [[ $denovo == "false" ]]; then
    printf '%s\n' "ERROR]: None of the methods, denovo or reference based, selected. >Quitting" >&2
    end_process
fi 

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "vsearch paired-end pipeline with multiple sequencing runs in multiRunDir"
    echo "Process = chimera filtering"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "qualFiltered_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = chimera filtering"
    DIRS=$(pwd)
    printf "\n workingDir = $DIRS \n"
fi

#############################
### Start of the workflow ###
#############################
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
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/chimeraFiltered_out"
        export output_dir

        ### Prepare working env and check paired-end data
        first_file_check
        prepare_SE_env
    else
        #output dir
        output_dir=$"/input/chimeraFiltered_out"
        export output_dir
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check paired-end data
        prepare_SE_env
    fi

    ### Process samples
    for file in *.$fileFormat; do
        ### Make temporary directory for temp files (for each sample)
        if [[ -d tempdir ]]; then
            rm -rf tempdir
        fi 
        mkdir -p tempdir
        #Read file name; without extension
        input=$(echo $file | sed -e "s/.$fileFormat//")
        ## Preparing files for the process
        printf "\n____________________________________\n"
        printf "Processing $input ...\n"
        #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        check_gz_zip_SE
        ### Check input formats (fastq/fasta supported)
        check_extension_fastx

        #If input is FASTQ then convert to FASTA
        if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
            checkerror=$(seqkit fq2fa -t dna --line-width 0 $input.$extension -o $input.fasta 2>&1)
            check_app_error
            printf "Note: converted $extension to FASTA \n"

            extension=$"fasta"
            export extension
            was_fastq=$"true"
            export was_fastq
        fi

        ###############################
        ### Start chimera filtering ###
        ###############################
        #make output dir for CHIMERAS
        mkdir -p $output_dir/chimeras
        #dereplicate sequences
        if [[ $denovo == "true" ]]; then
            checkerror=$(vsearch --derep_fulllength $input.fasta \
            $minuniquesize \
            --sizein --sizeout \
            --fasta_width 0 \
            --uc tempdir/$input.dereplicated.uc \
            --output tempdir/$input.derep.fasta 2>&1)
            check_app_error

            #pre-cluster sequences; sorts seqs automaticcaly by decreasing abundance
            checkerror=$(vsearch --cluster_size tempdir/$input.derep.fasta \
            $cores \
            $id \
            --strand both \
            --sizein --sizeout \
            --fasta_width 0 \
            --uc tempdir/$input.preclustered.uc \
            --centroids tempdir/$input.preclustered.fasta 2>&1)
            check_app_error

            #search chimeras
            checkerror=$(vsearch --uchime_denovo tempdir/$input.preclustered.fasta \
            $abskew \
            $minh \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --chimeras $output_dir/chimeras/$input.denovo.chimeras.fasta \
            --nonchimeras tempdir/$input.fasta 2>&1)
            check_app_error

            if [[ $reference_based == "undefined" ]]; then
                #Extract all non-chimeric sequences and add to $output_dir
                checkerror=$(vsearch --usearch_global $input.fasta \
                -db tempdir/$input.fasta \
                --sizein --xsize \
                $id \
                --strand both \
                --fasta_width 0 \
                --matched $output_dir/$input.fasta 2>&1)
                check_app_error

                #If input was fastq, then move all converted FASTA files to $output_dir/FASTA
                if [[ $was_fastq == "true" ]]; then
                    mkdir -p $output_dir/chimeraFilt_input_to_FASTA
                    mv $input.fasta $output_dir/chimeraFilt_input_to_FASTA
                fi
            else
                checkerror=$(vsearch --uchime_ref tempdir/$input.fasta \
                $cores \
                --db $database \
                --sizein \
                --sizeout \
                --fasta_width 0 \
                --chimeras $output_dir/chimeras/$input.ref.chimeras.fasta \
                --nonchimeras tempdir/$input.ref.denovo.nonchimeras.fasta 2>&1)
                check_app_error

                #Extract all non-chimeric sequences
                checkerror=$(vsearch --usearch_global $input.fasta \
                -db tempdir/$input.ref.denovo.nonchimeras.fasta \
                --sizein --xsize \
                $id \
                --strand both \
                --fasta_width 0 \
                --matched $output_dir/$input.fasta 2>&1)
                check_app_error

                #If input was fastq, then move all converted FASTA files to $output_dir/FASTA
                if [[ $was_fastq == "true" ]]; then
                    mkdir -p $output_dir/chimeraFilt_input_to_FASTA
                    mv $input.fasta $output_dir/chimeraFilt_input_to_FASTA
                fi
            fi
        
        else #only reference based chimera filtering
            checkerror=$(vsearch --uchime_ref $input.fasta \
            $cores \
            --db $database \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --chimeras $output_dir/chimeras/$input.ref.chimeras.fasta \
            --nonchimeras $output_dir/$input.fasta 2>&1)
            check_app_error

            #If input was fastq, then move all converted FASTA files to $output_dir/FASTA
            if [[ $was_fastq == "true" ]]; then
                mkdir -p $output_dir/chimeraFilt_input_to_FASTA
                mv $input.fasta $output_dir/chimeraFilt_input_to_FASTA
            fi
        fi
    done

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    printf "\nCleaning up and compiling final stats files ...\n"
    if [[ $was_fastq == "true" ]]; then
        #Delete tempdirs
        if [[ $debugger != "true" ]]; then
            if [[ -d tempdir ]]; then
                rm -rf tempdir
            fi
            if [[ -d tempdir2 ]]; then
                rm -rf tempdir2
            fi
        else 
            #compress files in /tempdir
            if [[ -d tempdir ]]; then
                pigz tempdir/*
            fi
        fi
        #make stats
        cd $output_dir/chimeraFilt_input_to_FASTA
        mkdir -p tempdir2
        clean_and_make_stats
        cd ..
        # no need for chimeraFilt_input_to_FASTA, remove it
        rm -rf cd $output_dir/chimeraFilt_input_to_FASTA
    else
        clean_and_make_stats
    fi
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file
    printf "# Chimera filtering with vsearch.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Files in 'chimeraFiltered_out' directory represent chimera filtered sequences.
Files in 'chimeraFiltered_out/chimeras' directory represent identified putative chimeric sequences.
If input was FASTQ formatted file(s), then it was converted to FASTA, and only FASTA is outputted.

    Core commands -> \n" > $output_dir/README.txt
    if [[ $denovo == "true" ]]; then
        printf "denovo filtering: vsearch --uchime_denovo input.preclustered.fasta $abskew $minh --chimeras chimeras/output.denovo.chimeras.fasta --nonchimeras output.fasta \n" >> $output_dir/README.txt
    fi
    if [[ $reference_based != "undefined" ]]; then
        printf "reference based filtering: vsearch --uchime_ref input.fasta $cores --db database_file --chimeras chimeras/output.ref.chimeras.fasta --nonchimeras output.fasta \n" >> $output_dir/README.txt
    fi

    printf "\nSummary of sequence counts in 'seq_count_summary.txt'

##############################################
###Third-party applications for this process:
#vsearch (version $vsearch_version) for chimera filtering
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
#seqkit (version $seqkit_version) for converting fastq to fasta (if input was fastq)
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
################################################" >> $output_dir/README.txt
    ### if working with multiRunDir then cd /input/multiRunDir
    if [[ $multiDir == "TRUE" ]]; then 
        cd /input/multiRunDir
    fi
done

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

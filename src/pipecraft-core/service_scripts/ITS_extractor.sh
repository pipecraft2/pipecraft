#!/bin/bash

#Input = single-end fasta/fastq files. FASTQ files will be converted to FASTA files; output is only FASTA.

# Extraction of ITS region(s)

################################################
###Third-party applications:
#ITSx v1.1.3
#seqkit v2.3.0
#mothur 1.46.1
#pigz v2.4
#perl v5.32.0
#python3 with biopython
################################################
printf "# pipeline = $pipeline\n"
printf "# service = $service for $organisms\n"

#load variables
organisms=$"-t ${organisms}"
regions=$"--save_regions ${regions}"
partial=$"--partial ${partial}"
cores=$"--cpu ${cores}"
eval=$"-E ${e_value}"
score=$"-S ${scores}"
domains=$"-N ${domains}" 
complement=${complement}
only_full=${only_full}
truncate=${truncate}
region_for_clustering=${region_for_clustering}
cluster_full_and_partial=${cluster_full_and_partial}

#additional options, if selection != undefined
if [[ $complement == "false" ]]; then
    complement_in=$"--complement F"
else
    complement_in=$"--complement T"
fi
if [[ $only_full == "false" ]]; then
    only_full_in=$"--only_full F"
else
    only_full_in=$"--only_full T"
fi
if [[ $truncate == "false" ]]; then
    truncate_in=$"--truncate F"
else
    truncate_in=$"--truncate T"
fi
if [[ $cluster_full_and_partial == "true" ]]; then
    full_and_partial=$"full_and_partial"
else
    full_and_partial=$""
fi


# Source for functions
source /scripts/submodules/framework.functions.sh
#python module for removing empty fasta records if using --partial
run_python_module=$"python3 /scripts/submodules/remove_empty_seqs.py"

# check if working with multiple runs or with a single sequencing run
if [[ -d "/input/multiRunDir" ]]; then
    echo "ITSx with multiple sequencing runs in multiRunDir"
    echo "Process = ITSx"
    echo "pipeline = $pipeline"
    cd /input/multiRunDir
    # read in directories (sequencing sets) to work with. Skip directories renamed as "skip_*"
    if [[ $pipeline == "vsearch_OTUs" ]]; then
        DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "chimeraFiltered_out" | grep -v "chimeraFiltered_out.dada2"| grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    fi
    if [[ $pipeline == "UNOISE_ASVs" ]]; then
        DIRS=$(find . -maxdepth 2 -mindepth 1 -type d | grep "qualFiltered_out" | grep -v "skip_" | grep -v "merged_runs" | grep -v "tempdir" | sed -e "s/^\.\///") 
    fi
    echo "working in dirs:"
    echo $DIRS
    multiDir=$"TRUE"
    export multiDir
else
    echo "Working with individual sequencing run"
    echo "Process = ITSx"
    echo "pipeline = $pipeline"
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
        output_dir=$"/input/multiRunDir/${seqrun%%/*}/ITSx_out"
        export output_dir

        ### Prepare working env and check single-end data
        first_file_check
        prepare_SE_env
    else
        #output dir
        output_dir=$"/input/ITSx_out"
        export output_dir
        # Check if files with specified extension exist in the dir
        first_file_check
        # Prepare working env and check single-end data
        prepare_SE_env
    fi

    #make output dir for no-detections
    mkdir $output_dir/no_detections
    ### Process samples
    for file in *.$fileFormat; do
        ### Make temporary directory for temp files (for each sample)
        if [ -d tempdir ]; then
            rm -rf tempdir
        fi
        mkdir tempdir
        #Read file name; without extension
        input=$(echo $file | sed -e "s/.$fileFormat//")
        ## Preparing files for the process
        printf "\n____________________________________\n"
        printf "Processing $input ...\n"
        #If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
        check_gz_zip_SE
        ### Check input formats (fastq supported)
        check_extension_fastx
        #If input is FASTQ then convert to FASTA
        if [[ $extension == "fastq" ]] || [[ $extension == "fq" ]]; then
            checkerror=$(seqkit fq2fa -t dna --line-width 0 $input.$extension -o $input.fasta 2>&1)
            check_app_error
            printf "Note: converted $extension to FASTA \n"

            extension=$"fasta"
            export extension
            was_fastq=$"TRUE"
        fi

        ##################
        ### Start ITSx ###
        ##################
        #dereplicate sequences
        checkerror=$(mothur "#unique.seqs(fasta=$input.$extension)" 2>&1)
        check_app_error
        mv $input.unique.$extension tempdir
        mv $input.names tempdir

        #Run ITSx
        echo "ITSx -i $input.$extension -o $input. --preserve T --graphical F $organisms  $partial $regions $cores $eval  $score $domains  $complement_in $only_full_in $truncate_in"

        checkerror=$(ITSx -i tempdir/$input.unique.$extension \
        -o tempdir/$input. \
        --preserve T \
        --graphical F \
        $organisms \
        $partial \
        $regions \
        $cores \
        $eval \
        $score \
        $domains \
        $complement_in \
        $only_full_in \
        $truncate_in 2>&1)
        check_app_error

        #Delete empty output files
        find tempdir -empty -type f -delete

        #map dereplicated sequences back to generate full output file
        if [[ -s tempdir/$input..SSU.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..SSU.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..SSU.redundant.fasta $output_dir
            mv $output_dir/$input..SSU.redundant.fasta $output_dir/$input.SSU.fasta
            mkdir -p $output_dir/SSU
            mv $output_dir/$input.SSU.fasta $output_dir/SSU
        fi
            #process for --partial SSU
        if [[ -s tempdir/$input..SSU.full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..SSU.full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..SSU.full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/SSU/full_and_partial
            $run_python_module $output_dir/$input..SSU.full_and_partial.redundant.fasta $output_dir/SSU/full_and_partial/$input.SSU.full_and_partial.fasta
            rm $output_dir/$input..SSU.full_and_partial.redundant.fasta
        fi

        if [[ -s tempdir/$input..ITS1.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..ITS1.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..ITS1.redundant.fasta $output_dir
            mv $output_dir/$input..ITS1.redundant.fasta $output_dir/$input.ITS1.fasta
            mkdir -p $output_dir/ITS1
            mv $output_dir/$input.ITS1.fasta $output_dir/ITS1
        fi
                #process for --partial ITS1
        if [[ -s tempdir/$input..ITS1.full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..ITS1.full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..ITS1.full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/ITS1/full_and_partial
            $run_python_module $output_dir/$input..ITS1.full_and_partial.redundant.fasta $output_dir/ITS1/full_and_partial/$input.ITS1.full_and_partial.fasta
            rm $output_dir/$input..ITS1.full_and_partial.redundant.fasta
        fi

        if [[ -s tempdir/$input..5_8S.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..5_8S.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..5_8S.redundant.fasta $output_dir
            mv $output_dir/$input..5_8S.redundant.fasta $output_dir/$input.5_8S.fasta
            mkdir -p $output_dir/5_8S
            mv $output_dir/$input.5_8S.fasta $output_dir/5_8S
        fi
                    #process for --partial 5_8S
        if [[ -s tempdir/$input..5_8S.full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..5_8S.full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..5_8S.full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/5_8S/full_and_partial
            $run_python_module $output_dir/$input..5_8S.full_and_partial.redundant.fasta $output_dir/5_8S/full_and_partial/$input.5_8S.full_and_partial.fasta
            rm $output_dir/$input..5_8S.full_and_partial.redundant.fasta
        fi

        if [[ -s tempdir/$input..ITS2.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..ITS2.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error      
            mv tempdir/$input..ITS2.redundant.fasta $output_dir        
            mv $output_dir/$input..ITS2.redundant.fasta $output_dir/$input.ITS2.fasta
            mkdir -p $output_dir/ITS2
            mv $output_dir/$input.ITS2.fasta $output_dir/ITS2
        fi
            #process for --partial ITS2
        if [[ -s tempdir/$input..ITS2.full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..ITS2.full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..ITS2.full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/ITS2/full_and_partial
            $run_python_module $output_dir/$input..ITS2.full_and_partial.redundant.fasta $output_dir/ITS2/full_and_partial/$input.ITS2.full_and_partial.fasta
            rm $output_dir/$input..ITS2.full_and_partial.redundant.fasta
        fi

        if [[ -s tempdir/$input..LSU.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..LSU.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..LSU.redundant.fasta $output_dir
            mv $output_dir/$input..LSU.redundant.fasta $output_dir/$input.LSU.fasta
            mkdir -p $output_dir/LSU
            mv $output_dir/$input.LSU.fasta $output_dir/LSU
        fi
                #process for --partial LSU
        if [[ -s tempdir/$input..LSU.full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..LSU.full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..LSU.full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/LSU/full_and_partial
            $run_python_module $output_dir/$input..LSU.full_and_partial.redundant.fasta $output_dir/LSU/full_and_partial/$input.LSU.full_and_partial.fasta
            rm $output_dir/$input..LSU.full_and_partial.redundant.fasta
        fi

        if [[ -s tempdir/$input..full.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..full.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..full.redundant.fasta $output_dir
            mv $output_dir/$input..full.redundant.fasta $output_dir/$input.full.fasta
            mkdir -p $output_dir/full_ITS
            mv $output_dir/$input.full.fasta $output_dir/full_ITS
        fi
            #process for --partial full ITS
        if [[ -s tempdir/$input..full_and_partial.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input..full_and_partial.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input..full_and_partial.redundant.fasta $output_dir
            mkdir -p $output_dir/full_ITS/full_and_partial
            $run_python_module $output_dir/$input..full_and_partial.redundant.fasta $output_dir/full_ITS/full_and_partial/$input.full_and_partial.fasta
            rm $output_dir/$input..full_and_partial.redundant.fasta
        fi
            #no detections
        if [[ -s tempdir/$input._no_detections.fasta ]]; then
            checkerror=$(mothur "#deunique.seqs(fasta=tempdir/$input._no_detections.fasta, name=tempdir/$input.names)" 2>&1)
            check_app_error
            mv tempdir/$input._no_detections.redundant.fasta $output_dir
            mv $output_dir/$input._no_detections.redundant.fasta $output_dir/no_detections/$input.no_detections.fasta
        fi
    done

    #################################################
    ### COMPILE FINAL STATISTICS AND README FILES ###
    #################################################
    printf "\nCleaning up and compiling final stats files ...\n"

    #for each output separately (SSU,ITS1,5_8S,ITS2,LSU,full)
        #file identifier string after the process
    if [[ -d $output_dir/SSU ]]; then
        outfile_addition=$"SSU"
        subdir=$"SSU"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/SSU/full_and_partial ]]; then
            outfile_addition=$"SSU.full_and_partial"
            subdir=$"SSU/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/ITS1 ]]; then
        outfile_addition=$"ITS1"
        subdir=$"ITS1"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/ITS1/full_and_partial ]]; then
            outfile_addition=$"ITS1.full_and_partial"
            subdir=$"ITS1/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/5_8S ]]; then
        outfile_addition=$"5_8S"
        subdir=$"5_8S"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/5_8S/full_and_partial ]]; then
            outfile_addition=$"5_8S.full_and_partial"
            subdir=$"5_8S/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/ITS2 ]]; then
        outfile_addition=$"ITS2"
        subdir=$"ITS2"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/ITS2/full_and_partial ]]; then
            outfile_addition=$"ITS2.full_and_partial"
            subdir=$"ITS2/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/LSU ]]; then
        outfile_addition=$"LSU"
        subdir=$"LSU"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/LSU/full_and_partial ]]; then
            outfile_addition=$"LSU.full_and_partial"
            subdir=$"LSU/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/full_ITS ]]; then
        outfile_addition=$"full"
        subdir=$"full_ITS"
        clean_and_make_stats_multidir
        if [[ -d $output_dir/full_ITS/full_and_partial ]]; then
            outfile_addition=$"full_ITS.full_and_partial"
            subdir=$"full_ITS/full_and_partial"
            subdir=$(echo $subdir | sed -e "s/\//\\\\\//g")
            clean_and_make_stats_multidir
        fi
    fi
    if [[ -d $output_dir/no_detections ]]; then
        outfile_addition=$"no_detections"
        subdir=$"no_detections"
        clean_and_make_stats_multidir
    fi

    #If initial input was FASTQ then remove converted FASTA files
    if [[ $was_fastq == "TRUE" ]]; then
        mkdir -p $output_dir/ITSx_input_to_FASTA
        mv *.fasta $output_dir/ITSx_input_to_FASTA
    fi

    if [[ $debugger != "true" ]]; then
        if [[ -d "tempdir2" ]]; then
            rm -rf tempdir2
        fi
    fi
    end=$(date +%s)
    runtime=$((end-start))

    #Make README.txt file
    printf "# ITS regions extracted with ITSx (see 'Core command' below for the used settings).

    Files in 'ITSx_out' directory represent sequences that passed ITS Extractor.
    Regions are placed under corrseponding directory (i.e., ITS2 sequences are in 'ITS2' directory).
    Files in /no_detections directory represent sequences where no ITS regions were identified.

    If input was FASTQ formatted file(s), then it was converted to FASTA, and only FASTA is outputted.
    Input FASTA files (converted from FASTQ) are in ITSx_out/ITSx_input_to_FASTA directory.

    Core command -> 
    ITSx -i input.unique.seqs -o output --preserve T --graphical F $organisms $partial $regions $cores $eval $score $domains $complement_in $only_full_in $truncate_in

    Total run time was $runtime sec.

    ##############################################
    ###Third-party applications for this process:
    #ITSx v1.1.3 for extracting ITS regions
        #citation: Bengtsson-Palme J., et al., 2013. Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods in Ecology and Evolution 4, 914-919.
        #microbiology.se/software/itsx/
    #seqkit v2.3.0 for converting fastq to fasta (if input was fastq)
        #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
        #https://bioinf.shenwei.me/seqkit/
    #mothur 1.46.1 for unique and deunique sequences prior and after extraction
        #citation: Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41
        #https://github.com/mothur/mothur
    ##############################################" > $output_dir/README.txt
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
if [[ $multiDir == "TRUE" ]]; then
    workingDir=$"/input/multiRunDir"
    echo "workingDir=$workingDir"
    # var for multiRunDir pipe
    printf "ITSx" > $workingDir/.prev_step.temp
    # ITSx_out subdirectory var_file for clustering for multiDir
    if [[ $region_for_clustering != "" ]]; then
        printf "$region_for_clustering" > $workingDir/.ITS_region_for_clustering.temp
        printf "$full_and_partial" > $workingDir/.ITS_full_and_partial.temp
    fi
else
    # ITSx_out subdirectory for clustering for individual seqrun
    if [[ $region_for_clustering != "" ]]; then
        echo "workingDir=$output_dir/$region_for_clustering/$full_and_partial"
    else
        echo "workingDir=$output_dir"
    fi
fi
echo "fileFormat=$extension"
echo "readType=single_end"






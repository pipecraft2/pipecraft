#!/bin/bash

# MetaWorks EVS pipeline for paired-end data

##########################################################
###Third-party applications:
#MetaWorks v1.11.1
    #citation: Porter, T.M., Hajibabaei, M. 2020. METAWORKS: A flexible, scalable bioinformatic pipeline for multi-marker biodiversity assessments. BioRxiv, doi: https://doi.org/10.1101/2020.07.14.202960.
    #Distributed under the GNU General Public License v3.0
    #https://github.com/terrimporter/MetaWorks
##########################################################

### TODO
# has to be gz files - if not then make to gz! unizip and make gz if needed or just make gz.
# test if metaworks cut primers reorients the reads! 
# validate primers  ok - -> also for cut primers!
################################################################

#Source for functions
#exec $SHELL
eval "$(conda shell.bash hook)"
source /scripts/submodules/framework.functions.sh
#output dir
output_dir="/input/metaworks_out"
mkdir -p $output_dir


#samples
extension=$fileFormat && export fileFormat  # must be gz files -> check with pigz at first !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
filename_structure=${filename_structure}
R1=$(echo $filename_structure | sed 's/R{read}/R1/')
R2=$(echo $filename_structure | sed 's/R{read}/R2/')


conda activate MetaWorks_v1.11.2

#ITSpart=$"ITS2" #or ITS1

#database for RDP
regex='[^/]*$'
db1_temp=$(echo $database | grep -oP "$regex")
db=$(printf "/extraFiles/$db1_temp")

#merge
quality_cutoff=${quality_cutoff}
min_overlap=${min_overlap}
mismatch_fraction=${mismatch_fraction}
match_fraction=${match_fraction}


#cut primers and quality filtering
fwd_tempprimer=${forward_primers}
rev_tempprimer=${reverse_primers}
minlen=${minlen}
qual_cutoff_5end=${qual_cutoff_5end}
qual_cutoff_3end=${qual_cutoff_3end}
maxNs=${maxNs}
error_rate=${error_rate}
#novaseq=${}
primer_overlap=${primer_overlap}

#printf "quality_cutoff=${quality_cutoff}
#min_overlap=${min_overlap}
#mismatch_fraction=${mismatch_fraction}
#match_fraction=${match_fraction}
#
#
##cut primers and quality filtering
#fwd_tempprimer=${forward_primers}
#rev_tempprimer=${reverse_primers}
#minlen=${minlen}
#qual_cutoff_5end=${qual_cutoff_5end}
#qual_cutoff_3end=${qual_cutoff_3end}
#maxNs=${maxNs}
#error_rate=${error_rate}
##novaseq=${}
#primer_overlap=${primer_overlap}"

#novaseq=$"yes" #when data is from NextSeq or Nova
# if [[ $novaseq == "yes" ]]; then
#     novaseq=$"--nextseq-trim=$qual_cutoff_5end"
# else
#     novaseq=$""
# fi


#Denoise
minsize=${minsize}
marker=${marker} 

#pseudogene filtering
pseudogene_filtering=$"yes" #yes/no list

cores=$"2"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check () {
count=$(ls -1 *.$extension 2>/dev/null | wc -l)
if [[ $count != 0 ]]; then 
    :
else
    printf '%s\n' "ERROR]: cannot find files with specified extension '$extension'
Please check the extension of your files and specify again.
>Quitting" >&2
    end_process
fi 
}
first_file_check
### Prepare working env and check paired-end data
prepare_PE_env () {
#Remove 'old' output_dir if exist and make new empty one
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir
#Make tempdir2, for seq count statistics
if [[ -d tempdir2 ]]; then
    rm -rf tempdir2
fi
mkdir -p tempdir2
#Make a file where to read R1 and R2 file names for paired-end read processing.
touch tempdir2/files_in_folder.txt
for file in *.$extension; do
    echo $file >> tempdir2/files_in_folder.txt
done
#Check for empty spaces in the files names. Replace with _
while read file; do
    if [[ $file == *" "* ]]; then
        printf '%s\n' "WARNING]: File $file name contains empty spaces, replaced 'space' with '_'" >&2
        rename 's/ /_/g' "$file"
    fi
done < tempdir2/files_in_folder.txt
#Fix also names in files_in_folder file if they contained space
sed -i 's/ /_/g' tempdir2/files_in_folder.txt
#Check if R1 string is in the file name (if so, then assume that then reverse file has R2 in the file name)
grep "R1" < tempdir2/files_in_folder.txt > tempdir2/paired_end_files.txt || true
    #Check if everything is ok considering file names
if [[ -s tempdir2/paired_end_files.txt ]]; then
    :
else
    printf '%s\n' "ERROR]: no paired-end read files found.
File names must contain 'R1' and 'R2' strings! (e.g. s01_R1.fastq and s01_R2.fastq)
>Quitting" >&2
    end_process
fi
#Check multiple occurrences of R1 and R2 strings (e.g. R123.R1.fq). 
while read file; do
    x=$(echo $file | grep -o -E '(R1|R2)' | wc -l)
    if [[ $x == "1" ]]; then
        :
    elif [[ $x == "0" ]]; then
        printf '%s\n' "ERROR]: $file name does not contain R1 or R2 strings to identify paired-end reads. Remove file from folder or fix the name.
>Quitting" >&2
        end_process
    else    
        printf '%s\n' "ERROR]: $file name contains multiple R1 or R2 strings -> change names (e.g. R123.R1.fastq to S123.R1.fastq)
>Quitting" >&2
        end_process
    fi
done < tempdir2/files_in_folder.txt
}
prepare_PE_env

### Make primers.fasta
fwd_primer_array=$(echo $fwd_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
rev_primer_array=$(echo $rev_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
i=1
for fwd_primer in $fwd_primer_array; do
    #validate primer seq
    printf ">primer\n$fwd_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
    check_app_error

    for rev_primer in $rev_primer_array; do
        #validate primer seq
        printf ">primer\n$rev_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
        check_app_error

        echo ">primer_pair$i" >> /input/tempdir2/primers.fasta
        #reverse complement rev primer
        rev_primer_rc=$(printf ">rev_primer\n$rev_primer" | seqkit seq --reverse --complement --seq-type DNA -v | sed -n 2p)
        echo "$fwd_primer...$rev_primer_rc" >> /input/tempdir2/primers.fasta
        ((i=i+1))
    done
done
primers=$"/input/tempdir2/primers.fasta"

## ENTER THE MetaWorks
cd /MetaWorks1.11.2

if [[ -f config_ESV.pipecraft.yaml ]]; then
    rm config_ESV.pipecraft.yaml
fi

# Make configuration file for MetaWorks v1.11.2
printf "# Configuration file for MetaWorks v1.11.2

# Author: Teresita M. Porter
# Date: August 30, 2022
# Slightly modified for PipeCraft2, by Sten Anslan 03.02.2023

############################################################################
# How to use MetaWorks v1
############################################################################

# 1. CREATE the MetaWorks_v1 environment, run this step one time only
# conda env create -f environment.yml

# 2. ACTIVATE the MetaWorks_v1 environment, always do this before running the pipeline
# conda activate MetaWorks_v1.11.2

# 3. EDIT the config.yaml file so customize directory names, paths, and choose your options for your MetaWorks run (below)

# 4. RUN snakemake
# snakemake --jobs 8 --snakefile snakefile_ESV --configfile config_ESV.yaml

# 5. DEACTIVATE the conda environment after the run is finished
# conda deactivate

# For instructions on installing conda, ORFfinder, or custom reference sets please refer to the MetaWorks website at https://terrimporter.github.io/MetaWorksSite/tutorial

############################################################################
# Identify raw read files

# This directory should contain the compressed paired-end Illumina reads, ex. *.fastq.gz
raw: "/input"

# Indicate 'sample' and 'read' wildcards from the raw filenames in the data folder (above):
# Ex. Sample filename structure:
# 	SITE-CONDITION-REPLICATE_S1_L001_R1_001.fastq.gz
# 	{sample}_L001_R{read}_001.fastq.gz
raw_sample_read_wildcards: "/input/{sample}_R{read}.fq.gz"

# SEQPREP sample wildcard and parameters
# These files should be in a data folder (above)
# Ex.
#	{sample}_L001_R1_001.fastq.gz
raw_sample_forward_wildcard: "/input/{sample}_R1.fq.gz"
raw_sample_reverse_wildcard: "/input/{sample}_R2.fq.gz"

############################################################################
# Directory for the output

# This directory will be created to contain pipeline results for a marker
# ex. COI, ITS, SSU
# keep the name short and simple with no spaces or weird punctuation, underscores are okay
dir: "/input/metaworks1.11.2_out"

############################################################################
# Raw read pairing

SEQPREP:
# Phred score quality cutoff (default 13):
    q: 13
# Minimum overlap (bp) length between forward and reverse reads:
    o: 25
# Maximum fraction of mismatches allowed in overlap (default 0.02):
    m: 0.02
# Minimum fraction of matching overlap (default 0.90):
    n: 0.90

############################################################################
# Primer trimming

# Each marker should have a fasta file with anchored linked adapters
# A single primer pair is sufficient or multiple sets if used for the same marker
# ex. COI_BE, COI_F230R, COI_mljg
# >AmpliconName;
# ^FwdPrimerSeq...ReverseComplementedRevPrimerSeq$
CUTADAPT:
    fasta: "/input/primers.fasta"

# Minimum sequence length (bp) to retain after trimming primers:
    m: 150
# Phred quality score cutoffs at the ends:
    q: "20,20"
# Error rate (default 0.1)
    e: 0.1
# Minimum adapter overlap (default 3)
    O: 24
# Maximum number of N's:
    mn: 0

############################################################################
# Denoising

# Indicate minimum number of reads per cluster to retain (default 8)
VSEARCH_DENOISE:
    minsize: 4

############################################################################
# ESV x sample table

VSEARCH_TABLE:
# Indicate number of threads to use
    t: 1

# Which marker classifier will you be using?
# Choose from ['16S', '18S_eukaryota', '18S_diatom', '12S_fish', '12S_vertebrate', 'ITS_fungi', '28S_fungi', 'rbcL_eukaryota', 'rbcL_diatom', 'rbcL_landPlant', 'ITS_plants', or 'COI']
marker: 'COI'

############################################################################
# ITSx extractor (edit if needed otherwise skip over this section)

# Indicate which spacer region to focus on:
# Choose from ['ITS1' or 'ITS2']
ITSpart: 'ITS2'

############################################################################
# Taxonomic assignment

RDP:
# enter the amount of memory to allocate to the RDP classifier here (default 8g):
    memory: "-Xmx100g"

# Do you want to use a custom-trained dataset?
# Set to 'yes' if using the following classifiers:
# COI, 12S_fish, 12S_vertebrate, rbcL_eukaryota, rbcl_diatom, 18S_eukaryota, 18S_diatom, ITS_UNITE, ITS_plants
# Set to 'no' if using RDP built-in classifiers:
# 16S or ITS_fungi (lsu or warcup)
# Choose from ['yes' or 'no']
    custom: 'yes'

# If you are using a custom-trained reference set 
# enter the path to the trained RDP classifier rRNAClassifier.properties file here:
    t: "/extraFiles2/rRNAClassifier.properties"
# If you are using the 16S RDP classifier built-in reference set, the pipeline will use these params:
    c: 0
    f: "fixrank"

# Otherwise you are using one of the RDP classifier built-in fungal classifiers:
# Choose from: ['fungallsu', 'fungalits_unite', 'fungalits_warcup']
    g: 'fungallsu'

###########################################################################
# Pseudogene filtering

# Indicate if you want to filter out putative pseudogenes:
# Set to 'no' for an rRNA gene/spacer region) or if working with protein coding gene but don't want to screen out putative pseudogenes then skip over the rest of this section
# Set to 'yes' if working with a protein coding gene and you want to screen out putative pseudogenes
# Choose from: ['yes' or 'no']
pseudogene_filtering: 'yes'

# Grep is used to refine the output to a single broad taxonomic group according to expected primer specificity 
# MetaWorks uses NCBI taxonomy, whole ranks only, see https://www.ncbi.nlm.nih.gov/taxonomy

# Simple grep search example:
# Ex. To target vertebrates do
## grep '-e Chordata'
# run this command with genetic code 2 below

# Compound grep search example:
# Ex. To target invertebrates do
## grep '-e Metazoa' to target Metazoa
## grep '-v Chordata' to exclude vertebrates
# run this with genetic code 5 below

# Grep search type [1 or 2] 
# (1) simple grep search (only taxon1 filter will be used) or 
# (2) a compound grep search (taxon1 and taxon2 filters will be used
grep_type: 1

taxon1: '-e Metazoa'
taxon2: '-v Chordata'

# If pseudogene_filtering was set to 'yes' then select pseudogene removal_type here [1|2]
# There are two pseudogene filtering methods available:
# (1) removal of sequences with unusually short/long open reading frames (ex. COI, rbcL)
# (2) HMM profile analysis and removal of sequences with unusually low HMM scores (ex. COI arthropoda)

removal_type: 1

# If removal_type was set to 2, then indicate the name of the hmm profile (only available for COI arthropoda at this time)
hmm: 'bold.hmm'

# Translate ESVs into all open reading frames
ORFFINDER:

# genetic code:
# 1 = standard code (use for rbcL)
# 2 = vertebrate mitochondrial (use for COI if targeting vertebrates)
# 5 = invertebrate mitochondrial (use for COI if targeting invertebrates)
# See NCBI for additional genetic codes:
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    g: 5

# ORF start codon to use:
# 0 = ATG only
# 1 = ATG and alternative initiation codon (ORFfinder default)
# 2 = any sense codon
    s: 2

# minimum length (ORFfinder default 75, min 30 nt)
    ml: 300

# ignore nested ORFs [true or false]
    n: 'true'

# strand [both, plus, minus]
    strand: 'plus'

###########################################################################
# Output options

# Option 1: If you have tens to low hundreds of samples, 
# you can choose to print a single combined report:
# results.csv lists all ESVs per sample with read counts and assigned taxonomy in a single file

# Option 2: If you have high hundreds to thousands of samples, 
# it is more memory- and time-efficient to work with the component files separately
# The separate ESV.table, taxonomy.csv, and sequence/extracted ITS/ORF FASTA files are all indexed by the same ESV id (Zotu #) and can be combined later in R for integrated analyses.  These files will be listed in results.csv .

# report_type [1 or 2]
report_type: 1
" > /input/config_ESV.pipecraft.yaml


#########################################################################################################
#########################################################################################################
#########################################################################################################
# Make snakefile_ESV for pipecraft
cd /input/ 
if [[ -f snakefile_ESV.pipecraft ]]; then
    rm snakefile_ESV.pipecraft
fi

cp /MetaWorks1.11.2/snakefile_ESV snakefile_ESV.pipecraft
sed -i 's/Read in vars from config_ESV.yaml/Read in vars from config_ESV.pipecraft.yaml/' snakefile_ESV.pipecraft
sed -i '6 i\configfile: "config_ESV.pipecraft.yaml"' snakefile_ESV.pipecraft

cd /MetaWorks1.11.2
# Run MetaWorks
echo "Running snakemake"
#checkerror=$(
    snakemake --jobs 2 --snakefile /input/snakefile_ESV.pipecraft --configfile /input/config_ESV.pipecraft.yaml
    # 2>&1) 
#check_app_error

echo "Snakemake done"

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fq.gz"

echo "readType=paired_end"


# ###########################
# ### create an ESV table ###
# ###########################
# # The results.csv file is the final file in the MetaWorks pipeline. 
# # R using reshape2 library

# # Run R in MetaWorks output folder

# # Read the MetaWorks results
# results <- read.csv("results.csv", header = TRUE)
# # Reshape the results into an 'OTU table'
# ESVtable <- reshape2::dcast(results, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum) 
    # ## giver error if results.csv is empty 
     ### Error in dim(ordered) <- ns : dims [product 1] do not match the length of object [0]

# ##########################
# ### get final taxonomy ###
# ##########################
# # replace _ with " " in the results.csv data frame
# results_sub <- data.frame(lapply(results, function(x) {gsub("_", " ", x)}))
# #convert ESVsize back to integer
# results_sub$ESVsize = as.integer(results_sub$ESVsize)
# #reshape funn table with taxonomy
# full_table = reshape2::dcast(results_sub, SampleName ~ GlobalESV + ORFseq + Root + rBP + SuperKingdom + skBP + Kingdom + kBP + Phylum + pBP + Class + cBP + Order + oBP + Family + fBP + Genus + gBP + Species + sBP,  
#                                 value.var = "ESVsize", fun.aggregate = sum)
# # set ',' as a separator for OTU + Seq + Tax + bootstrap values
# names(full_table) = stringr::str_replace_all(names(full_table), "_", ",")
# #transpose full table
# t_full_table = t(full_table)
# #take only tax info and discard first row (=SampleName)
# tax = rownames(t_full_table)[2:nrow(as.data.frame(t_full_table))] 
# #convert to data.frame for tidyr
# df.tax = as.data.frame(tax, col.names = "tax")
# #split tax columns
# tax_table = tidyr::separate(df.tax, col = "tax", sep = ",", into = c("GlobalESV", "ORFseq", 
#                                                                "Root", "rBP", 
#                                                                "SuperKingdom", "skBP", 
#                                                                "Kingdom", "kBP", 
#                                                                "Phylum", "pBP", 
#                                                                "Class", "cBP", 
#                                                                "Order", "oBP", 
#                                                                "Family", "fBP", 
#                                                                "Genus", "gBP", 
#                                                                "Species", "sBP", NA), fill = "right")

# # Write tables to output and order records
# ESVtable_t = t(ESVtable)
# write.table(ESVtable_t[order(rownames(ESVtable_t)),],"final_ESV_table.txt", row.names = TRUE, col.names = FALSE)
# write.table(tax_table[order(tax_table$GlobalESV),],"final_tax_table.txt", row.names = FALSE, col.names = TRUE)

#write out ALSO ASVs.fasta!!! 

# #DONE
# # Final outputs = final_ESV_table.txt & final_tax_table.txt

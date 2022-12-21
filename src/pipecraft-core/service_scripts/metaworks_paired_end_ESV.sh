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
source /scripts/submodules/framework.functions.sh
#output dir
output_dir=$"/input/metaworks_out"

#samples
extension=$fileFormat # must be gz files -> check with pigz at first
filename_structure=${filename_structure}
R1=$(echo $filename_structure | sed 's/R{read}/R1/')
R2=$(echo $filename_structure | sed 's/R{read}/R2/')

echo "wd =" 
pwd

echo "activate CONDA"
exec $SHELL
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.11.2

#mkdir -p /.cache && chmod a+rwX /.cache
chmod 777 / 
#a=$(ls -la /opt/)
#printf "\n $a \n"

#ITSpart=$"ITS2" #or ITS1
#db=$"/files/CO1v4_trained/rRNAClassifier.properties" #path to RDP database -> make this as default in the container?


#merge
quality_cutoff=${quality_cutoff}
min_overlap=${min_overlap}
mismatch_fraction=${mismatch_fraction}
match_fraction=${match_fraction}
novaseq=$"yes" #when data is from NextSeq or Nova

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

if [[ $novaseq == "yes" ]]; then
    novaseq=$"--nextseq-trim=$qual_cutoff_5end"
else
    novaseq=$""
fi


#Denoise
minsize=${minsize}
marker=${marker} 

#pseudogene filtering
pseudogene_filtering=$"yes" #yes/no list

cores=$"6"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
#prepare_PE_env

# mkdir /input/tempdir2

# ### Make primers.fasta
# fwd_primer_array=$(echo $fwd_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
# rev_primer_array=$(echo $rev_tempprimer | sed 's/,/ /g' | sed 's/I/N/g')
# i=1
# for fwd_primer in $fwd_primer_array; do
#     #validate primer seq
#     printf ">primer\n$fwd_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
#     check_app_error

#     for rev_primer in $rev_primer_array; do
#         #validate primer seq
#         printf ">primer\n$rev_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
#         check_app_error

#         echo ">primer_pair$i" >> /input/tempdir2/primers.fasta
#         #reverse complement rev primer
#         rev_primer_rc=$(printf ">rev_primer\n$rev_primer" | seqkit seq --reverse --complement --seq-type DNA -v | sed -n 2p)
#         echo "$fwd_primer...$rev_primer_rc" >> /input/tempdir2/primers.fasta
#         ((i=i+1))
#     done
# done
# primers=$"/input/primers.fasta"

### ENTER THE MetaWorks
# cd /MetaWorks1.11.1

# if [[ -f config_ESV.pipecraft.yaml ]]; then
#     rm config_ESV.pipecraft.yaml
# fi

# # Make configuration file for MetaWorks v1.11.1
# printf "# Configuration file for MetaWorks v1.11.1

# # Author: Teresita M. Porter
# # Date: August 26, 2022
# # This is slightly modified version for PipeCraft2, by Sten Anslan (September 2022). Core processes of MetaWorks are not edited.
# ############################################################################
# # Identify raw read files

# # This directory should contain the compressed paired-end Illumina reads, ex. *.fastq.gz
# raw: '/input/'

# # Indicate 'sample' and 'read' wildcards from the raw filenames in the data folder (above):
# raw_sample_read_wildcards: '/input/$filename_structure.$extension'

# # SEQPREP sample wildcard and parameters
# raw_sample_forward_wildcard: '/input/$R1.$extension'
# raw_sample_reverse_wildcard: '/input/$R2.$extension'

# ############################################################################
# # Directory for the output

# # This directory will be created to contain pipeline results
# dir: 'input/metaworks_out'

# ############################################################################
# # Raw read pairing

# SEQPREP:
# # Phred score quality cutoff (default 13):
#     q: $quality_cutoff
# # Minimum overlap (bp) length between forward and reverse reads:
#     o: $min_overlap
# # Maximum fraction of mismatches allowed in overlap (default 0.02):
#     m: $mismatch_fraction
# # Minimum fraction of matching overlap (default 0.90):
#     n: $match_fraction

# ############################################################################
# # Primer trimming

# # Each marker should have a fasta file with anchored linked adapters
# # A single primer pair is sufficient or multiple sets if used for the same marker
# # ex. COI_BE, COI_F230R, COI_mljg
# # >AmpliconName;
# # ^FwdPrimerSeq...ReverseComplementedRevPrimerSeq$
# CUTADAPT:
#     fasta: '$primers'

# # Minimum sequence length (bp) to retain after trimming primers:
#     m: $minlen
# # Phred quality score cutoffs at the ends:
#     q: '$qual_cutoff_5end,$qual_cutoff_3end'
# # Error rate (default 0.1)
#     e: $error_rate
# # Minimum adapter overlap
#     O: $primer_overlap
# # Maximum number of N's:
#     mn: $maxNs

# ############################################################################
# # Denoising

# # Indicate minimum number of reads per cluster to retain
# # Default here is to remove all singletons and doubletons and retain clusters with 3+ reads
# VSEARCH_DENOISE:
#     minsize: $minsize

# ############################################################################
# # ESV x sample table

# VSEARCH_TABLE:
# # Indicate number of threads to use
#     t: 4

# # Which marker classifier will you be using?
# # Choose from ['16S', '18S_eukaryota', '18S_diatom', '12S_fish', '12S_vertebrate', 'ITS_fungi', '28S_fungi', 'rbcL_eukaryota', 'rbcL_diatom', 'rbcL_landPlant', 'ITS_plants', or 'COI']
# marker: '$marker'

# ############################################################################
# # ITSx extractor (edit if needed otherwise skip over this section)

# # Indicate which spacer region to focus on:
# # Choose from ['ITS1', 'ITS2']
# ITSpart: 'ITS2'

# ############################################################################
# # Taxonomic assignment

# RDP:
# # enter the amount of memory to allocate to the RDP classifier here (default 8g):
#     memory: '-Xmx8g'

# # Do you want to use a custom-trained dataset?
# # Set to 'yes' if using the following classifiers:
# # COI, 12S_fish, 12S_vertebrate, rbcL_eukaryota, rbcl_diatom, 18S_eukaryota, 18S_diatom, ITS_UNITE, ITS_plants
# # Set to 'no' if using RDP built-in classifiers:
# # 16S or ITS_fungi (lsu or warcup)
# # Choose from ['yes' or 'no']
#     custom: 'no'

# # If you are using a custom-trained reference set 
# # enter the path to the trained RDP classifier rRNAClassifier.properties file here:
#     t: ''

# # If you are using the 16S RDP classifier built-in reference set, the pipeline will use these params:
#     c: 0
#     f: 'fixrank'

# # Otherwise you are using one of the RDP classifier built-in fungal classifiers:
# # Choose from: ['fungallsu', 'fungalits_warcup']
#     g: 'fungallsu'

# ###########################################################################
# # Pseudogene filtering

# # Indicate if you want to filter out putative pseudogenes:
# # Set to 'no' for an rRNA gene/spacer region) or if working with protein coding gene but don't want to screen out putative pseudogenes then skip over the rest of this section
# # Set to 'yes' if working with a protein coding gene and you want to screen out putative pseudogenes
# # Choose from: ['yes' or 'no']
# pseudogene_filtering: 'no'

# # Grep is used to refine the output to a single broad taxonomic group according to expected primer specificity 
# # MetaWorks uses NCBI taxonomy, whole ranks only, see https://www.ncbi.nlm.nih.gov/taxonomy
# # To keep one taxon do '-e Chordata' with genetic code 2 to target vertebrates
# # To keep more than one taxon do '-e taxon1 -e taxon2' etc.
# # To exclude a taxon use '-v taxon3'
# # To keep one taxon and exlude a sub-taxon
# # ex. To target invertebrates with genetic code 5
# ## '-e Metazoa rdp.out.tmp | grep -v Chordata'
# taxon: '-e Arthropoda'

# # If pseudogene_filtering was set to 'yes' then select pseudogene removal_type here
# # There are two pseudogene filtering methods available:
# # (1) removal of sequences with unusually short/long open reading frames (ex. rbcL)
# # (2) HMM profile analysis and removal of sequences with unusually low HMM scores

# removal_type: 1

# # If removal_type was set to 2, then indicate the name of the hmm profile (only available for COI arthropoda at this time)
# hmm: 'bold.hmm'

# # Translate ESVs into all open reading frames
# ORFFINDER:

# # genetic code:
# # 1 = standard code (use for rbcL)
# # 2 = vertebrate mitochondrial (use for COI if targeting vertebrates)
# # 5 = invertebrate mitochondrial (use for COI if targeting invertebrates)
# # See NCBI for additional genetic codes:
# # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#     g: 5

# # ORF start codon to use:
# # 0 = ATG only
# # 1 = ATG and alternative initiation codon (ORFfinder default)
# # 2 = any sense codon
#     s: 2

# # minimum length (ORFfinder default 75, min 30 nt)
#     ml: 30

# # ignore nested ORFs (true|false)
#     n: 'true'

# # strand (both|plus|minus)
#     strand: 'plus'

# ###########################################################################
# # Output options

# # Option 1: If you have tens to low hundreds of samples, 
# # you can choose to print a single combined report:
# # results.csv lists all ESVs per sample with read counts and assigned taxonomy in a single file

# # Option 2: If you have high hundreds to thousands of samples, 
# # it is more memory- and time-efficient to work with the component files separately
# # The separate ESV.table, taxonomy.csv, and sequence/extracted ITS/ORF FASTA files are all indexed by the same ESV id (Zotu #) and can be combined later in R for integrated analyses.  These files will be listed in results.csv .

# # report_type (1|2)
# report_type: 1
# " > /input/tempdir2/config_ESV.pipecraft.yaml

# cd /input/ 
ls 
pwd
#########################################################################################################
#########################################################################################################
#########################################################################################################

# Make snakefile_ESV for pipecraft
# if [[ -f snakefile_ESV.pipecraft ]]; then
#     rm snakefile_ESV.pipecraft
# fi

#cp snakefile_ESV /input/tempdir2/snakefile_ESV.pipecraft
#sed -i 's/Read in vars from config_ESV.yaml/Read in vars from config_ESV.pipecraft.yaml/' /input/tempdir2/snakefile_ESV.pipecraft
#sed -i '6 i\configfile: "config_ESV.pipecraft.yaml"' /input/tempdir2/snakefile_ESV.pipecraft


# Run MetaWorks
echo "Running snakemake"
#checkerror=$(snakemake -h 2>&1)
#check_app_error

checkerror=$(snakemake --jobs 2 --snakefile snakefile_ESV --configfile config_ESV.pipecraft.yaml 2>&1) 
#/input/tempdir2/snakefile_ESV.pipecraft 
check_app_error

echo "Snakemake done"

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Summary of sequence counts in '$output_dir/seq_count_summary.txt'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=fq.gz"
echo "dataFormat=$dataFormat"
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

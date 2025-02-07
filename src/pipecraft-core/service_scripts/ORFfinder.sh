#!/bin/bash

#Using ORFfinder to identify open reading frames + arthropod COI HMM based filtering
#Workflow based on implementations in MetaWorks v1.12.0 (by Teresita M. Porter)
# Input = fasta file.

##################################
###Third-party applications:
#ORFfinder
#seqkit
#pigz
##################################
# input fasta file
    regex='[^/]*$'
    rep_seqs_file_path=$(echo $fasta_file | grep -oP "$regex")
    rep_seqs_temp=$(basename $rep_seqs_file_path) #basename, needed for macOS
    rep_seqs_file=$(printf "/extraFiles/$rep_seqs_temp")
    echo "rep_seqs_file = $rep_seqs_file"

# min_length    = positive integer. minimum length of the output sequence
# max_length    = positive integer. maximum length of the output sequence
# genetic_code  = positive integer. genetic code for translation. 5 = invertebrate mitochondrial code. Specify values from 1 to 33
# start_codon   = list [0, 1 or 2]. 0 = ATG only; 1 = ATG and alternative initation codons; 2 = any sense codon
# ignore_nested = bool. TRUE = ignore nested ORFs (completely placed within another)
# strand        = list [plus, minus, both]. output ORFs on specified strand only

#output dir
output_dir=$"/input/"
# Source for functions
source /scripts/submodules/framework.functions.sh

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.12.0

### Make ORFfinder tempdir
if [[ -d "$output_dir/tempdir_filter_numts" ]]; then
    rm -rf $output_dir/tempdir_filter_numts
fi
mkdir -p $output_dir/tempdir_filter_numts

# Check if repfasta_file_seqs is fasta
extension=$(echo $fasta_file | (awk 'BEGIN{FS=OFS="."} {print $NF}';))
if [[ $extension == "fasta" ]] || [[ $extension == "fa" ]] || [[ $extension == "fas" ]]; then
    :
else
    printf '%s\n' "ERROR]: $extension formatting not supported here for 'rep seqs file'!
Supported extensions: fasta, fas, fa.
>Quitting" >&2
    end_process
fi

### Using ORFfinder to identify open reading frames
echo "# Running ORFfinder | "
    # -outfmt 0 = list of ORFs in FASTA format (aa)
    # -outfmt 1 = CDS in FASTA format (nucl)
checkerror=$(ORFfinder -in $rep_seqs_file \
    -ml $min_length \
    -g $genetic_code \
    -s $start_codon \
    -n $ignore_nested \
    -strand $strand \
    -outfmt 1 \
    -out $output_dir/tempdir_filter_numts/ORFs.nt.fasta 2>&1)
check_app_error

### retain the longest ORF
echo "# Sorting ORFfinder -outfmt 1 results | "
in_name=$(basename $rep_seqs_file | awk -F "." '(NF=NF-1)' | sed -e 's/ /\./')

#below parsing does not work if sequence IDs have _
perl /MetaWorks1.12.0/perl_scripts/parse_orfs3.plx $output_dir/tempdir_filter_numts/ORFs.nt.fasta $output_dir/$in_name.ORFs.temp

# filter by max_length
checkerror=$(seqkit seq --quiet -g -w 0  -M $max_length $output_dir/$in_name.ORFs.temp \
    > $output_dir/$in_name.ORFs.fasta 2>&1)
check_app_error

# remove temp file
rm $output_dir/$in_name.ORFs.temp

# make list of ORFs
checkerror=$(seqkit seq --quiet --name $output_dir/$in_name.ORFs.fasta \
                            > $output_dir/$in_name.ORFs.list.txt 2>&1)
check_app_error

# write putative pseudogenes to a file
checkerror=$(seqkit grep --quiet -w 0 -v -f $output_dir/$in_name.ORFs.list.txt $rep_seqs_file \
                                                    > $output_dir/$in_name.notORFs.fasta 2>&1)
check_app_error

# make list of notORFs
checkerror=$(seqkit seq --quiet -n $output_dir/$in_name.notORFs.fasta \
                                    > $output_dir/$in_name.notORFs.list.txt 2>&1)
check_app_error

echo "# ORFfinder DONE | "

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
if [[ $debugger != "true" ]]; then
    rm -r $output_dir/tempdir_filter_numts
fi
# count outputs
input_seqs=$(grep -c "^>" $rep_seqs_file)
nonORFs=$(grep -c "^>" $output_dir/$in_name.notORFs.fasta)
ORFs=$(grep -c "^>" $output_dir/$in_name.ORFs.fasta)

end=$(date +%s)
runtime=$((end-start))

#Make README.txt file for demultiplexed reads
in=$(echo $in_name | sed -e 's/\/extraFiles\///')
printf "# Used ORFfinder to remove putative pseudogenes and off-targets. 
    (details in the MetaWorks user guide: https://terrimporter.github.io/MetaWorksSite/details/ - sequences are translated into every possible open reading frame (ORF) using ORFfinder, the longest ORF is reatined. Putative pseudogenes are removed as outliers with unusually small/large ORF lengths. Outliers are calcualted as follows: Sequence lengths shorter than the 25th percentile - 1.5*IQR (inter quartile range) are removed as putative pseudogenes (or sequences with errors that cause a frame-shift). Sequence lengths longer than the 75th percentile + 1.5*IQR are also removed as putative pseudogenes.

Input file $in.$extension contained $input_seqs sequences.
Filtered output file $in.ORFs.fasta contains $ORFs sequences.

Generated files:
# $in.ORFs.fasta       = filtered output by the ORFfinder; sequences with the longest ORFs (open reading frames); contains $ORFs reads.
# $in.ORFs.list.txt    = list of sequence identifiers of the above file.
# $in.notORFs.fasta    = putative pseudogenes/off-target ORFs; contains $nonORFs reads.
# $in.notORFs.list.txt = list of sequence identifiers of the above file.

Core commands -> 
ORFfinder -in $rep_seqs_file_path -ml $min_length -g $genetic_code -s $start_codon -n $ignore_nested -strand $strand
seqkit seq --quiet -g -w 0  -M $max_length $in_name.ORFs.temp > $in_name.ORFs.fasta \n" > $output_dir/README_ORFfinder.txt

printf "Total run time was $runtime sec.

##############################################
###Third-party applications for this process:
# MetaWorks v1.12.0 (strategy for filtering putative NUMTs via ORFfinder)
    citation: Porter, T. M., & Hajibabaei, M. (2022). MetaWorks: A flexible, scalable bioinformatic pipeline for high-throughput multi-marker biodiversity assessments. PLoS One, 17(9), e0274260. https://doi.org/10.1371/journ al.pone.0274260
# ORFfinder v0.4.3 for finding ORFs
    #https://www.ncbi.nlm.nih.gov/orffinder/
##############################################" >> $output_dir/README_ORFfinder.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"
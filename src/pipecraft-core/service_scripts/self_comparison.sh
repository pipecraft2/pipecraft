#!/bin/bash

# Make pairwise comparison for sequences in a fasta file
 # method = vsearch/blast
 # input = fasta file
 # output = pairwise comparison file

################################################
###Third-party applications:
# BLAST 2.16.0+
# vsearch v2.29.4
##############################################

# Checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
blast_version=$(blastn -version 2>&1 | head -n 1 | awk '{print $2}')
printf "# Checking tool versions ...\n"
printf "# vsearch (version $vsearch_version)\n"
printf "# BLAST (version $blast_version)\n"


#get specified input fasta file
regex='[^/]*$'
input_fasta_temp=$(echo $fasta_file | grep -oP "$regex")
input_fasta=$(printf "/extraFiles/$input_fasta_temp")
printf "\n input fasta = $input_fasta \n"

#Source for functions
source /scripts/submodules/framework.functions.sh

#output dir
output_dir=$"/input"

#############################
### Start of the workflow ###
#############################
#start time
start_time=$(date)
start=$(date +%s)

# fasta basename
fasta_basename=$(basename $input_fasta)

# Check if input is a fasta file
printf "# Checking if input is in fasta format...\n"
if ! head -n1 "$input_fasta" | grep -q "^>"; then
    log_error "Input file does not appear to be in fasta format (should start with >)"
    exit 1
fi



### Pairwise comparison
# Pairwise comparison with BLAST
if [[ $method == "BLAST" ]]; then
    printf "\n#Making blast database from the input fasta \n"
    checkerror=$(makeblastdb -in $input_fasta -parse_seqids -dbtype nucl 2>&1)
    check_app_error

    printf "# Pairwise comparison with BLASTn \n"
    checkerror=$(blastn -db $input_fasta \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand' \
            -out $output_dir/$fasta_basename.self_comp.blast.temp \
            -qcov_hsp_perc $coverage \
            -perc_identity $identity \
            -query $input_fasta \
            -num_threads $cores \
            -strand $strand 2>&1)
            check_app_error

    # add header to BLAST output
    printf "# Adding header to BLAST output\n"
    sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tqcovhsp\tsstrand' $output_dir/$fasta_basename.self_comp.blast.temp

    # Remove self-hits from BLAST output
    printf "# Removing self-hits from BLAST output\n" 
    checkerror=$(awk '$1 != $2' $output_dir/$fasta_basename.self_comp.blast.temp > $output_dir/$fasta_basename.self_comp.blast.txt 2>&1)
    check_app_error

    # Remove temporary BLAST output file
    rm $output_dir/$fasta_basename.self_comp.blast.temp
fi


# Pairwise comparison with vsearch
if [[ $method == "vsearch" ]]; then
    printf "# Pairwise comparison with vsearch --usearch_global \n"
    #convert perc_identity and coverage_perc for vsearch
    vsearch_perc_identity=$(awk "BEGIN {print $identity/100}")
    vsearch_coverage_perc=$(awk "BEGIN {print $coverage/100}")
    checkerror=$(vsearch --usearch_global $input_fasta \
            --db $input_fasta \
            --strand $strand \
            --self \
            --id $vsearch_perc_identity \
            --iddef 2 \
            --userout $output_dir/$fasta_basename.self_comp.vsearch.txt \
            --userfields query+target+id+alnlen+qcov+tcov+ql+tl+ids+mism+gaps+qilo+qihi+qstrand+tstrand \
            --maxaccepts 0 \
            --maxrejects 0 \
            --query_cov $vsearch_coverage_perc \
            --threads $cores 2>&1)
    check_app_error

    # add header to vsearch output
    printf "# Adding header to vsearch output\n"
    sed -i '1i query\ttarget\tid\talnlen\tqcov\ttcov\tql\ttl\tids\tmism\tgaps\tqilo\tqihi\tqstrand\ttstrand' $output_dir/$fasta_basename.self_comp.vsearch.txt
fi

########################################
### CLEAN UP AND COMPILE README FILE ###
########################################
for i in $input_fasta.n*; do
    [[ -f $i ]] || continue
    rm -f "$i"
done

end=$(date +%s)
runtime=$((end-start))

###Make README.txt file
if [[ $method == "BLAST" ]]; then
    pairwise_comarison=$"blastn -db $fasta_basename -outfmt 6 -out $fasta_basename.self_comp.blast.temp -qcov_hsp_perc $coverage -perc_identity $identity -query $input_fasta -num_threads $cores -strand $strand

qseqid   - Query sequence ID
sseqid   - Subject sequence ID
pident   - Percentage identity
length   - Alignment length
mismatch - Number of mismatches
gapopen  - Number of gap openings
qstart   - Start of alignment in query
qend     - End of alignment in query
sstart   - Start of alignment in subject
send     - End of alignment in subject
evalue   - Expect value
bitscore - Bit score
qlen     - Query sequence length
slen     - Subject sequence length
qcovs    - Query coverage per subject
qcovhsp  - Query coverage per HSP
sstrand  - Subject strand"
fi

if [[ $method == "vsearch" ]]; then
    pairwise_comarison=$"vsearch --usearch_global $fasta_basename --db $fasta_basename --strand $strand --self --id $vsearch_perc_identity --iddef 2 --userout $fasta_basename.self_comp.vsearch.txt --userfields query+target+id+alnlen+qcov+tcov+ql+tl+ids+mism+gaps+qilo+qihi+qstrand+tstrand --maxaccepts 0 --maxrejects 0 --query_cov $vsearch_coverage_perc --threads $cores

query   - Query sequence ID
target  - Subject sequence ID
id      - Percentage identity (according to --iddef 2)
alnlen  - Alignment length (set to 0 if there is no alignment)
qcov    - Percentage of the query sequence that is aligned with the target sequence. Computed as 100.0 * (matches + mismatches) / query sequence length. Internal or terminal gaps are not taken into account.
tcov    - Percentage of the target sequence that is aligned with the query sequence. Computed as 100.0 * (matches + mismatches) / target sequence length. Internal or terminal gaps are not taken into account.
ql      - Query sequence length
tl      - Target sequence length
ids     - Number of matches in the alignment
mism    - Number of mismatches in the alignment
gaps    - Number of gaps in the alignment
qilo    - First nucleotide of the query aligned with the target (ignoring initial gaps). Nucleotide numbering starts from 1.
qihi    - Last nucleotide of the query aligned with the target (ignoring terminal gaps). Nucleotide numbering starts from 1.
qstrand - Query strand orientation
tstrand - Subject strand orientation
"
fi

printf "# Performed pairwise comparison with $method (see 'Core commands' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Core commands -> 
pairwise comparison: $pairwise_comarison

################################################
###Third-party applications:
# BLAST (version $blast_version) (if BLAST was used)
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. 
# vsearch (version $vsearch_version) (if vsearch was used)
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
##############################################" > $output_dir/README.self_comparison.txt

#Done
printf "\nDONE "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$fileFormat"
echo "readType=$readType"


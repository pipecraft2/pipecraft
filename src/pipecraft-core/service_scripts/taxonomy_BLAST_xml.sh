#!/bin/bash

#Taxonomy annotation with BLAST
#Input = single-end fasta file + database file.
#If database = fasta file, then makeblastdb with BLAST+, otherwise run BLAST
#Outputs: BLAST_1st_best_hit.txt, BLAST_10_best_hits.txt

##########################################################
###Third-party applications:
    #BLAST 2.12.0+
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421. 
#python3 with biopython module
#gawk
##########################################################

#load variables
extension=$fileFormat

#database [# here, a possibility for multiple databases to be added: database=$"-db $db1 $db2 $db3 $db4 $db5"]
regex='[^/]*$'
db1_temp=$(echo $database_file | grep -oP "$regex")
db1=$(printf "/extraFiles/$db1_temp")
echo "db1 = $db1"
#mandatory options
task=$"-task ${task}" # list: blastn, megablast
strands=$"-strand ${strands}" #list: both, plus
#additional options
cores=$"-num_threads ${cores}" # positive integer
evalue=$"-evalue=${e_value}" # float
wordsize=$"-word_size=${word_size}" # positive integer
reward=$"-reward=${reward}" # positive integer
penalty=$"-penalty=${penalty}" # negative integer
gapopen=$"-gapopen=${gap_open}" # positive integer
#gapextend setting (default=undefined with megablast)
if [[ $gap_extend == null ]] || [[ -z $gap_extend ]] || [[ $gap_extend == "undefined" ]]; then
    gapextend=$""
else
    gapextend=$"-gapextend=${gap_extend}" 
fi

# Source for functions
source /scripts/submodules/framework.functions.sh
parseXML=$"/scripts/submodules/Blast_xml_parse.py"
#output dir
output_dir=$"/input/taxonomy_out"

#############################
### Start of the workflow ###
#############################
start=$(date +%s)
### Check if files with specified extension exist in the dir
first_file_check
### Prepare working env and check paired-end data
prepare_SE_env
#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
	
check_gz_zip_SE
### Check input formats (fasta supported)
check_extension_fasta
### Select last fasta file in the folder as input for BLAST
for file in *.$extension; do
	IN=$(echo $file)
done
echo "input = $IN"

### Check and assign BLAST database
d1=$(echo $db1 | awk 'BEGIN{FS=OFS="."}{print $NF}') #get the extension
#make blast database if db is not formatted for BLAST
db_dir=$(dirname $db1)
check_db_presence=$(ls -1 $db_dir/*.nhr 2>/dev/null | wc -l)
if [[ $check_db_presence != 0 ]]; then
	if [[ $d1 == "fasta" ]] || [[ $d1 == "fa" ]] || [[ $d1 == "fas" ]] || [[ $d1 == "fna" ]] || [[ $d1 == "ffn" ]]; then
		database=$"-db $db1"
	elif [[ $d1 == "ndb" ]] || [[ $d1 == "nhr" ]] || [[ $d1 == "nin" ]] || [[ $d1 == "not" ]] || [[ $d1 == "nsq" ]] || [[ $d1 == "ntf" ]] || [[ $d1 == "nto" ]]; then
		db1=$(echo $db1 | awk 'BEGIN{FS=OFS="."}NF{NF-=1};1')
		database=$"-db $db1"
	fi
elif [[ $d1 == "fasta" ]] || [[ $d1 == "fa" ]] || [[ $d1 == "fas" ]] || [[ $d1 == "fna" ]] || [[ $d1 == "ffn" ]]; then
		printf '%s\n' "Note: converting fasta formatted database for BLAST"
		makeblastdb -in $db1 -input_type fasta -dbtype nucl
		database=$"-db $db1"
fi

## Perform taxonomy annotation
printf '%s\n' "Running BLAST"
checkerror=$(blastn \
-query $IN \
$strands \
$cores \
$database \
-out $output_dir/10BestHits.xml \
$task \
-max_target_seqs 10 \
$evalue \
$wordsize \
$reward \
$penalty \
$gapopen \
$gapextend \
-max_hsps 1 \
-outfmt=5 2>&1)
check_app_error

### Parse BLAST xml file
cd $output_dir
python3 /scripts/submodules/Blast_xml_parse.py 10BestHits.xml

#Extract first best hit from 10 best hits file
gawk 'BEGIN{FS=OFS="\t"}{print $1,$4}' < 10hits.txt | \
sed -e 's/^$/NO_BLAST_HIT/g' | \
sed '1 i\SeqID+first_hit+score+e-value+query len+query start+query end+target len+target start+target end+align len+identities+gaps+coverage%+id%' \
> BLAST_1st_best_hit.txt

#Format 10 best hits from 10 best hits file
gawk 'BEGIN{FS=OFS="\t"}{print $1,$4,$6,$7,$9,$10,$12,$13,$15,$16,$18,$19,$21,$22,$24,$25,$27,$28,$30,$31}' < 10hits.txt | \
sed '1 i\SeqID+first_hit+score+e-value+query len+query start+query end+target len+target start+target end+align len+identities+gaps+coverage%+id%+next_hits' \
> BLAST_10_best_hits.txt

cd ..

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"
rm $output_dir/10hits.txt
rm $output_dir/10BestHits.xml
if [[ -d tempdir2 ]];then
	rm -r tempdir2
fi

end=$(date +%s)
runtime=$((end-start))


#Make README.txt file
printf "Taxonomy annotation with BLAST. 
Input = $IN
BLAST_1st_best_hit.txt = BLAST results for the 1st best hit in the used database(s).
BLAST_10_best_hits.txt = BLAST results for the 10 best hits in the used database(s).

BLAST values filed separator is '+'. When pasting the taxonomy results to e.g. Excel, then first denote '+' as as filed separator to align the columns.

score -> blast score
e-value -> blast e-value (note that e-value of 0 in BLAST means e-value is < 5e-324).
query len -> query (i.e. OTU/ASV) sequence length
query start -> start position of match in the query seq
query end -> end position of match in the query seq
target len -> target seq length in the database
target start -> start position of match in the target seq
target end -> end position of match in the target seq
align len -> alignment length of query and target
identities -> number of identical matches
gaps -> number of gaps in the alignment
coverage -> query coverage percentage against the target sequence (100 percent is full-length match; low coverage may indicate presence of chimeric sequence/OTU)
id -> identity percentage against the target sequence

Total run time was $runtime sec.

##########################################################
###Third-party applications [PLEASE CITE]:
    #BLAST 2.11.0+
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. 
#####################################################" > $output_dir/README.txt

#Done
printf "\nDONE\n"
printf "Data in directory '$output_dir'\n"
printf "Check README.txt files in output directory for further information about the process.\n"
printf "Total time: $runtime sec.\n\n"

#variables for all services
echo "workingDir=$output_dir"
echo "fileFormat=$extension"

echo "readType=single_end"

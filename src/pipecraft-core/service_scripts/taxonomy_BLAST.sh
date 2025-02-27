#!/bin/bash

#Taxonomy annotation with BLAST
#Input = single-end fasta file + database file.
#If database = fasta file, then makeblastdb with BLAST+, otherwise run BLAST
#Outputs: BLAST_1st_best_hit.txt, BLAST_10_best_hits.txt

##########################################################
###Third-party applications:
# BLAST+
# python3 with biopython module
# seqkit
# gawk
##########################################################

#load variables
regex='[^/]*$'
db1_temp=$(echo $database_file | grep -oP "$regex")
db1=$(printf "/extraFiles/$db1_temp")
echo "db1 = $db1"

### input fasta file
IN=$(echo $fasta_file | grep -oP "$regex")
# overwrite fileFormat variable; get it from input fasta_file
fileFormat=$(echo $IN | awk -F. '{print $NF}')
export fileFormat

#mandatory options
task=$"-task ${task}"       # e.g. blastn, megablast
strands=$"-strand ${strands}"   # both, plus

#additional options
cores=$"-num_threads ${cores}"  # positive integer
evalue=$"-evalue=${e_value}"    # float
wordsize=$"-word_size=${word_size}" # positive integer
reward=$"-reward=${reward}"     # positive integer
penalty=$"-penalty=${penalty}"  # negative integer
gapopen=$"-gapopen=${gap_open}" # positive integer

#gapextend setting (default=undefined with megablast)
if [[ $gap_extend == null ]] || [[ -z $gap_extend ]] || [[ $gap_extend == "undefined" ]]; then
    gapextend=$""
else
    gapextend=$"-gapextend=${gap_extend}"
fi

# Source for functions
source /scripts/submodules/framework.functions.sh

#output dir
output_dir=$"/input/taxonomy_out"

#############################
### Start of the workflow ###
#############################
echo "output_dir = $output_dir"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir

### Start time
start_time=$(date)
start=$(date +%s)

### Check if files with specified extension exist in the dir
first_file_check

#If input is compressed, then decompress (keeping the compressed file, but overwriting if filename exists!)
check_gz_zip_SE

### Check input formats (fasta supported)
check_extension_fasta

### Select last fasta file in the folder as input for BLAST
for file in *.$fileFormat; do
    IN=$(echo $file)
done
echo "input = $IN"

### Check and assign BLAST database
d1=$(echo $db1 | awk 'BEGIN{FS=OFS="."}{print $NF}') #get the extension
db_dir=$(dirname $db1)
check_db_presence=$(ls -1 $db_dir/*.nhr 2>/dev/null | wc -l)

# If the database is not yet formatted for BLAST, format it if it's in FASTA form
if (( $check_db_presence != 0 )); then
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
    printf '%s\n' "BLAST database was converted to BLAST format (file = $db1)."
fi

## Do BLAST
printf '%s\n' "Running BLAST"
checkerror=$(blastn \
-query $IN \
$strands \
$cores \
$database \
-out $output_dir/10BestHits.txt \
$task \
-max_target_seqs 10 \
$evalue \
$wordsize \
$reward \
$penalty \
$gapopen \
$gapextend \
-max_hsps 1 \
-outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" 2>&1)
check_app_error

### Parse 10BestHits ### QUERY SEQUENCE HEADERS MUST BE UNIQUE!
cd $output_dir

#get only first occurrence of a duplicate row (1st hit)
mkdir -p tempdir
awk 'BEGIN{FS="+"} !seen[$1]++' 10BestHits.txt > tempdir/1.temphit
#check which seqs got a hit
gawk 'BEGIN{FS="+"}{print $1}' < tempdir/1.temphit | uniq > tempdir/gothits.names

#add no_hits flag to sequences that didn't get a hit
seqkit seq -n ../$IN > tempdir/$IN.names
grep -v -w -F -f tempdir/gothits.names tempdir/$IN.names | sed -e 's/$/\tNo_significant_similarity_found/' >> tempdir/1.temphit

#add header
sed -e '1 i\qseqid+1st_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident' tempdir/1.temphit > tempdir/BLAST_1st_hit.txt

### For hits 2-10
for i in {2..10}; do
    awk -v i="$i" 'BEGIN{FS="+"} ++seen[$1]==i' 10BestHits.txt > tempdir/$i.temphit
    gawk 'BEGIN{FS="+"}{print $1}' < tempdir/$i.temphit | uniq > tempdir/gothits.names
    grep -v -w -F -f tempdir/gothits.names tempdir/$IN.names | sed -e 's/$/\tNo_BLAST_hit/' >> tempdir/$i.temphit && rm tempdir/gothits.names
done

# Sort each .temphit and merge
for file in tempdir/*.temphit; do
    sort -k 1 --field-separator=+ "$file" > "$file.temp" && rm "$file"
done

paste \
  tempdir/1.temphit.temp \
  tempdir/2.temphit.temp \
  tempdir/3.temphit.temp \
  tempdir/4.temphit.temp \
  tempdir/5.temphit.temp \
  tempdir/6.temphit.temp \
  tempdir/7.temphit.temp \
  tempdir/8.temphit.temp \
  tempdir/9.temphit.temp \
  tempdir/10.temphit.temp > BLAST_10_best_hits.txt

rm tempdir/*.temp

# Clean up the final 10-best-hits file
sed -i 's/No_significant_similarity_found.*/No_significant_similarity_found/' BLAST_10_best_hits.txt
sed -i 's/No_BLAST_hit.*/No_BLAST_hit/' BLAST_10_best_hits.txt
sed -i '1 i\qseqid+1st_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+2nd_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+3rd_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+4th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+5th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+6th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+7th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+8th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+9th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+10th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident' BLAST_10_best_hits.txt

##### BLAST 1st hit with query SEQ ######
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ../$IN | sed -e 's/\r//' > tempdir/$IN.oneline
sort -k 1 --field-separator=\t tempdir/BLAST_1st_hit.txt > BLAST_1st_best_hit.temp
sed -i 's/qseqid.*//' BLAST_1st_best_hit.temp
sed -i '/^$/d' BLAST_1st_best_hit.temp
sort -k 1 --field-separator=\t tempdir/$IN.oneline | sed -e 's/^>//' | sed -e 's/\r//' > tempdir/seqs.txt

paste tempdir/seqs.txt BLAST_1st_best_hit.temp > BLAST_1st_best_hit.txt && rm BLAST_1st_best_hit.temp
sed -i '1 i\qseqid+query_seq+qseqid+1st_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident' BLAST_1st_best_hit.txt

mv 10BestHits.txt tempdir/
sed -i 's/\t/+/g' BLAST_1st_best_hit.txt
sed -i 's/\t/+/g' BLAST_10_best_hits.txt

##########################################################
### ADD sim_score PER BLOCK IMMEDIATELY AFTER pident   ###
##########################################################
# 1) For BLAST_1st_best_hit.txt, we simply append sim_score at the end (after pident).
awk -F'+' 'BEGIN{OFS="+"}
NR==1 {
   print $0,"sim_score";
   next
}
{
   # pident is $19, length is $12, qlen is $5
   if($5+0 > 0){
     sim = $19 * ($12 / $5)
   } else {
     sim = "NA"
   }
   print $0, sim
}' BLAST_1st_best_hit.txt > BLAST_1st_best_hit.tmp && mv BLAST_1st_best_hit.tmp BLAST_1st_best_hit.txt

# 2) For BLAST_10_best_hits.txt, each block of 17 columns becomes 18 by inserting sim_score after pident.
awk -F'+' 'BEGIN{OFS="+"}
NR==1 {
   # Rebuild header: each block is 17 columns, we add sim_score after pident
   split($0, head, FS)
   nblocks = length(head)/17
   new_header = ""
   for(i=0; i<nblocks; i++){
     start = i*17 + 1
     end   = start + 16
     for(j=start; j<=end; j++){
       new_header = new_header head[j] OFS
     }
     # Insert "sim_score" label
     new_header = new_header "sim_score"
     if(i < nblocks-1){
       new_header = new_header OFS
     }
   }
   print new_header
   next
}
{
   split($0, a, FS)
   nblocks = length(a)/17
   row_out = ""
   for(i=0; i<nblocks; i++){
     start = i*17 + 1
     end   = start + 16
     # parse columns for sim_score
     qlen      = a[start+2]
     align_len = a[start+9]
     pident    = a[start+16]
     if(qlen+0>0){
       sim = pident * (align_len / qlen)
     } else {
       sim = "NA"
     }
     block_str = ""
     for(j=start; j<=end; j++){
       block_str = block_str a[j] OFS
     }
     # Insert sim_score
     block_str = block_str sim
     if(i < nblocks-1){
       block_str = block_str OFS
     }
     row_out = row_out block_str
   }
   print row_out
}' BLAST_10_best_hits.txt > BLAST_10_best_hits.tmp && mv BLAST_10_best_hits.tmp BLAST_10_best_hits.txt

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ...\n"

if [[ $debugger != "true" ]]; then
    if [[ -d tempdir ]]; then
        rm -r tempdir
    fi
    if [[ -d tempdir2 ]]; then
        rm -r tempdir2
    fi
fi

end=$(date +%s)
runtime=$((end-start))

db_x=$(echo $db1 | sed -e 's/\/extraFiles\///')

# Make README.txt file
printf "# Taxonomy was assigned using BLAST (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Query    = $IN
Database = $db_x

BLAST_1st_best_hit.txt = BLAST results for the 1st best hit in the used database (sim_score appended).
BLAST_10_best_hits.txt = BLAST results for the 10 best hits in the used database, with sim_score inserted after pident in each block.

BLAST values field separator is '+'. When pasting the taxonomy results to e.g. Excel, then first denote '+' as the field separator to align the columns.

qseqid    = query id
query_seq = query sequence
1st_hit   = first BLAST hit
qlen      = query sequence length
slen      = subject sequence length
qstart    = start of alignment in query
qend      = end of alignment in query
sstart    = start of alignment in subject
send      = end of alignment in subject
evalue    = expect value
length    = alignment length
nident    = number of identical matches
mismatch  = number of mismatches
gapopen   = number of gap openings
gaps      = total number of gaps
sstrand   = subject strand
qcovs     = query coverage per subject
pident    = percentage of identical matches
sim_score = similarity score of a hit taking the query coverage into account; calculated as (pident * (alignment length / qlen))

Core command ->
blastn -query $IN $strands $db_x $task -max_target_seqs 10 $evalue $wordsize $reward $penalty $gapopen $gapextend -max_hsps 1

Total run time was $runtime sec.

##########################################################
###Third-party applications [PLEASE CITE]:
    #BLAST 2.14.0+
    #citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421.
#####################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

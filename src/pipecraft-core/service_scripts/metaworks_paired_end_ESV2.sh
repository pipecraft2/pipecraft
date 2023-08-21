#!/bin/bash

### TODO
# DONE! has to be gz files - if not then make to gz! unizip and make gz if needed or just make gz.
# OK! test if metaworks cut primers reorients the reads! YES, but only in v1.12
# validate primers  ok - -> also for cut primers!

# Kui 16S, siis kasutab built-in DB, muidu peab alla laadima.. TESTi 16S (ITS warcupi ei kasuta, AINT UNITE)
# 16S Pro vs mt16S for euk. Test marker=?

#genetic_code DISABLED when pseudogene_filt = FALSE
#ITS region = DISABLED when other markers selected than ITS

# Check if filename structure files are found in the folder. 
# If trimmed size = 0 then stop right away.
################################################################
eval "$(conda shell.bash hook)"
conda activate MetaWorks_v1.12.0

set -e 

#Source for functions
source /scripts/submodules/framework.functions.sh

cd /MetaWorks1.12.0

extension=$fileFormat # must be gz files
filename_structure=${filename_structure}
R1=$(echo $filename_structure | sed "s/R{read}/R1/")
R2=$(echo $filename_structure | sed "s/R{read}/R2/")
primers=$"/input/primers.in.fasta" #compiled below

regex='[^/]*$'
db=$(echo $database | grep -oP "$regex")
database=$(basename $db)

marker=${marker}
quality_cutoff=${quality_cutoff}
min_seq_len=${min_seq_len}
minsize=${minsize}
min_overlap=${min_overlap}
match_fraction=${match_fraction}
mismatch_fraction=${mismatch_fraction}
qual_cutoff_3end=${qual_cutoff_3end}
qual_cutoff_5end=${qual_cutoff_5end}
primer_mismatch=${primer_mismatch}
primer_overlap=${primer_overlap}
maxNs=${maxNs}
pseudogene_filtering=${pseudogene_filtering}
genetic_code=${genetic_code}
orf_len=${orf_len}


if [[ $pseudogene_filtering == "TRUE" ]] || [[ $pseudogene_filtering == "true" ]]; then
    pseudogene_filtering=$"yes"
else
    pseudogene_filtering=$"no"
fi

### Check if input fq files are compressed. If not, then gz compress.
check_compress=$(echo $extension | (awk 'BEGIN{FS=OFS="."} {print $NF}';))
    if [[ $check_compress == "gz" ]]; then
        extension=$(echo $extension | (awk 'BEGIN{FS=OFS="."} {print $(NF-1), $NF}';))
        export extension
    elif [[ $check_compress == "fastq" ]]; then
        printf '%s\n' "WARNING]: Compressing fastq files for MetaWorks"
        pigz *.fastq
        extension=$"fastq.gz"
        export extension
    elif [[ $check_compress == "fq" ]]; then
        printf '%s\n' "WARNING]: Compressing fq files for MetaWorks"
        pigz *.fq
        extension=$"fq.gz"
        export extension
    else 
        printf '%s\n' "ERROR]: Please input fastq, fq (or .gz comressed) files for MetaWorks.
        >Quitting" >&2 
        quit_process  
    fi 

### Make primers fasta file
if [[ -f /input/primers.in.fasta ]]; then
    rm /input/primers.in.fasta
fi 
fwd_primer_array=$(echo $forward_primers | sed 's/,/ /g' | sed 's/I/N/g')
rev_primer_array=$(echo $reverse_primers | sed 's/,/ /g' | sed 's/I/N/g')
i=1
for fwd_primer in $fwd_primer_array; do
    #validate primer seq
    printf ">primer\n$fwd_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
    check_app_error

    for rev_primer in $rev_primer_array; do
        #validate primer seq
        printf ">primer\n$rev_primer" | checkerror=$(seqkit seq --seq-type DNA -v 2>&1)
        check_app_error

        echo ">primer_pair$i" >> /input/primers.in.fasta
        #reverse complement rev primer
        rev_primer_rc=$(printf ">rev_primer\n$rev_primer" | seqkit seq --reverse --complement --seq-type DNA -v | sed -n 2p)
        echo "$fwd_primer...$rev_primer_rc" >> /input/primers.in.fasta
        ((i=i+1))
    done
done


#####################################
### Edit variables in config file ###
#####################################
#input dir (should contain the compressed paired-end Illumina reads, ex. *.fastq.gz)
sed -e 's/raw: "testing\/COI_data"/raw: "\/input"/' config_ESV.yaml > in.config_ESV.yaml

#input file structure (filename_structure)
sed -i "s/raw_sample_read_wildcards:.*/raw_sample_read_wildcards: \"\/input\/$filename_structure.$extension\"/" in.config_ESV.yaml
sed -i "s/raw_sample_forward_wildcard:.*/raw_sample_forward_wildcard: \"\/input\/$R1.$extension\"/" in.config_ESV.yaml
sed -i "s/raw_sample_reverse_wildcard:.*/raw_sample_reverse_wildcard: \"\/input\/$R2.$extension\"/" in.config_ESV.yaml
#output dir (keep the name short and simple with no spaces or weird punctuation, underscores are okay)
sed -i 's/dir:.*/dir: "\/input\/metaworks_out"/' in.config_ESV.yaml
# Phred score quality cutoff (default 19):
sed -i "s/q: 13/q: $quality_cutoff/" in.config_ESV.yaml
# Minimum overlap (bp) length between forward and reverse reads (default 25)
sed -i "s/o: 25/o: $min_overlap/" in.config_ESV.yaml
# Maximum fraction of mismatches allowed in overlap (default 0.02)
sed -i "s/m: 0.02/m: $match_fraction/" in.config_ESV.yaml
# Minimum fraction of matching overlap (default 0.90)
sed -i "s/n: 0.90/n: $mismatch_fraction/" in.config_ESV.yaml
#primers file
sed -i 's/fasta: "testing\/adapters_anchored.fasta"/fasta: "\/input\/primers.in.fasta"/' in.config_ESV.yaml

# Minimum sequence length (bp) to retain after trimming primers
sed -i "s/m: 150/m: $min_seq_len/" in.config_ESV.yaml
# Phred quality score cutoffs at the ends:
sed -i "s/q: \"20,20\"/q: \"$qual_cutoff_3end,$qual_cutoff_5end\"/" in.config_ESV.yaml
# Error rate (default 1)
sed -i "s/e: 0.1/e: $primer_mismatch/" in.config_ESV.yaml
# Minimum adapter overlap (default 15)
sed -i "s/O: 3/O: $primer_overlap/" in.config_ESV.yaml
# Maximum number of N's
sed -i "s/mn: 3/mn: $maxNs/" in.config_ESV.yaml
# Enable reverse complement option [Yes, No] (default "Yes") - not adjustable in PipeCraft2 GUI
sed -i 's/rc: "No"/rc: "Yes"/' in.config_ESV.yaml
# Indicate minimum number of reads per cluster to retain (default 8)
sed -i "s/minsize: 8/minsize: $minsize/" in.config_ESV.yaml
# Indicate number of threads to use (vsearch)
sed -i "s/t: 15/t: 1/" in.config_ESV.yaml
# Which marker classifier will you be using? Choose from ['16S', '18S_eukaryota', '18S_diatom', '12S_fish', '12S_vertebrate', 'ITS_fungi', '28S_fungi', 'rbcL_eukaryota', 'rbcL_diatom', 'rbcL_landPlant', 'ITS_plants', or 'COI']
sed -i "s/marker: 'COI'/marker: '$marker'/" in.config_ESV.yaml
# ITSx extractor; choose from ['ITS1' or 'ITS2']
sed -i "s/ITSpart: 'ITS2'/ITSpart: '$ITS_region'/" in.config_ESV.yaml

### RDP
# amount of memory to allocate to the RDP classifier here (replace default 8g with 16g):
sed -i 's/memory: "-Xmx8g"/memory: "-Xmx16g"/' in.config_ESV.yaml
# Do you want to use a custom-trained dataset? Choose from ['yes' or 'no'] - "no" only for 16S
if [[ $marker == "16S" ]]; then
    sed -i "s/custom: 'yes'/custom: 'no'/" in.config_ESV.yaml
fi
# specify database (database) if using a custom-trained reference set 
sed -i "s/t: \"\/path\/to\/rRNAClassifier.properties\"/t: \"\/extraFiles2\/$database\"/" in.config_ESV.yaml

### Pseudogene filtering
sed -i "s/pseudogene_filtering: 'yes'/pseudogene_filtering: '$pseudogene_filtering'/" in.config_ESV.yaml
# minimum length of ORF
sed -i "s/ml: 30/ml: $orf_len/" in.config_ESV.yaml


#copy config file to WD
cp in.config_ESV.yaml /input



#run snakemake
snakemake --jobs $cores --snakefile snakefile_ESV --configfile in.config_ESV.yaml

echo "workingDir=$workingDir"
echo "fileFormat=$fileFormat"

echo "readType=paired_end"

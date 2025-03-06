#!/bin/bash

function convert_IUPAC () {
    local primer=$1
    local IUPAC_codes="R|Y|S|W|K|M|B|D|H|V|N|I"

    if echo $primer | grep -q -E "$IUPAC_codes"; then
        # define IUPAC codes and their corresponding sets of bases
        declare -A IUPAC_bases=(
            ["R"]="AG"
            ["Y"]="CT"
            ["S"]="GC"
            ["W"]="AT"
            ["K"]="GT"
            ["M"]="AC"
            ["B"]="CGT"
            ["D"]="AGT"
            ["H"]="ACT"
            ["V"]="ACG"
            ["N"]="ATCG"
            ["I"]="ATCG"
        )

        # replace IUPAC codes with random bases
        converted_primer=$(echo $primer | sed -e "s/\[//g; s/\]//g")
        for i in $(seq 1 ${#converted_primer}); do
            if echo ${converted_primer:i-1:1} | grep -q -E "$IUPAC_codes"; then
                IUPAC_code=${converted_primer:i-1:1}
                IUPAC_bases_set=${IUPAC_bases[$IUPAC_code]}
                random_base=${IUPAC_bases_set:$RANDOM%${#IUPAC_bases_set}:1}
                converted_primer=${converted_primer:0:i-1}${random_base}${converted_primer:i}
            fi
        done

        # return converted primer
        echo $converted_primer
    else
        # return original primer when no IUPAC codes were detected
        echo $primer
    fi
}

fwd_primer_array=$(echo $1 | sed 's/,/ /g' | sed 's/I/N/g')
rev_primer_array=$(echo $2 | sed 's/,/ /g' | sed 's/I/N/g')

# Remove existing output files
rm -f fwd_primer.fasta rev_primer.fasta fwd_primer_RC.fasta rev_primer_RC.fasta liked_fwd_revRC.fasta liked_rev_fwdRC.fasta

# Forward primer(s) to fasta file
i=1
for primer in $fwd_primer_array; do
    converted_primer=$(convert_IUPAC $primer)
    echo ">fwd_primer$i" >> fwd_primer.fasta
    echo $converted_primer >> fwd_primer.fasta
    ((i=i+1))
done

# Reverse primer(s) to fasta file
i=1
for primer in $rev_primer_array; do
converted_primer=$(convert_IUPAC $primer)
echo ">rev_primer$i" >> rev_primer.fasta
echo $converted_primer >> rev_primer.fasta
((i=i+1))
done

checkerror=$(seqkit seq --quiet -t dna -r -p fwd_primer.fasta >> fwd_primer_RC.fasta 2>&1)
checkerror=$(seqkit seq --quiet -t dna -r -p rev_primer.fasta >> rev_primer_RC.fasta 2>&1)

i=1
while read LINE; do
fwd_primer=$(echo $LINE | grep -v ">")
if [ -z "$fwd_primer" ]; then
:
else
while read LINE; do
rev_primer_RC=$(echo $LINE | grep -v ">")
if [ -z "$rev_primer_RC" ]; then
:
else
echo ">primer$i" >> liked_fwd_revRC.fasta
echo "$fwd_primer...$rev_primer_RC" >> liked_fwd_revRC.fasta
((i=i+1))
fi
done < rev_primer_RC.fasta
fi
done < fwd_primer.fasta

i=1
while read LINE; do
rev_primer=$(echo $LINE | grep -v ">")
if [ -z "$rev_primer" ]; then
:
else
while read LINE; do
fwd_primer_RC=$(echo $LINE | grep -v ">")
if [ -z "$fwd_primer_RC" ]; then
:
else
echo ">primer$i" >> liked_rev_fwdRC.fasta
echo "$rev_primer...$fwd_primer_RC" >> liked_rev_fwdRC.fasta
((i=i+1))
fi
done < fwd_primer_RC.fasta
fi
done < rev_primer.fasta




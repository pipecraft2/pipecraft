#!/bin/bash

## Script to compress homopolymer stretches (e.g., prior to k-mer counting)
## All homopolymer stretches will be collapsed to max H length (default, 1)

## Input:
# $1 = Parameter H (max homopolymer length)
# $2 = input FASTA file

## Usage example:
# 
# cat > input.fasta <<'EOT'
# >s1
# AACCCTTTTGGGGG
# >s2
# ACCTTTGGGGGGGG
# >s3
# ACTGACTGACTGAC
# EOT
# 
# homopolymer_compression.sh 2 input.fasta



awk -v H="$1" '\

BEGIN {
    if (H < 1) H = 1;
}

# If the line is a header, print it as is
/^>/ {
    print;
    next;
}

# Process sequence lines
{
    sequence = $0;
    collapsedSeq = "";
    count = 1;

    for (i = 2; i <= length(sequence); i++) {
        if (substr(sequence, i, 1) == substr(sequence, i - 1, 1)) {
            count++;
        } else {
            collapsedSeq = collapsedSeq substr(sequence, i - count, (count > H) ? H : count);
            count = 1;
        }
    }

    # Handle the last homopolymer stretch
    collapsedSeq = collapsedSeq substr(sequence, length(sequence) - count + 1, (count > H) ? H : count);

    print collapsedSeq;
}
' "$2"


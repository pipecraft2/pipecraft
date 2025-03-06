#!/usr/bin/env python

#Remove records in fasta file with empty sequences.
#ITSx reports empty sequences when using --partial option, use this script to clean the file.
#usage: script.py input.fasta output.fasta

from Bio import SeqIO
from sys import argv

print("python module running")
script, fasta_file, output = argv

with open(output, "w") as output_handle:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.seq != '':
            SeqIO.write(record, output_handle, "fasta-2line")
            
print("python module finished")
#convert fasta file to fastq file with dummy (40) quality scores

from Bio import SeqIO

with open("F3D2.merged_ASVs.rerep.fastq", "w") as output_handle:
    for record in SeqIO.parse("F3D2.merged_ASVs.rerep.fasta", "fasta"):
        record.letter_annotations["solexa_quality"] = [40] * len(record)
        SeqIO.write(record, output_handle, "fastq")

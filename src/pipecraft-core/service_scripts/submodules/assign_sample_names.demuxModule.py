#!/usr/bin/python

#Renaming files when using combinatorial dual indexes for demultiplexing, using cutadapt

from Bio import SeqIO
from sys import argv
import os

print("python module running")

#usage: script.py $indexes_file
script, indexes_file = argv

#read file extension from bash variable
extension = os.environ["fileFormat"]

#read demultiplexed file names
demux_files = open("tempdir2/demux_files.txt", 'r')
for line in demux_files:
	el = line.strip()
	el = el.split("-")

	#el[0] forward index name
	#el[1] reverse index name

	#open fwd index file and seach index seq for matching index name
	with open("tempdir2/index_fwd.uniq.renamed.fasta") as indexes_fwd:
		for fwdrecord in SeqIO.parse(indexes_fwd, "fasta"):
			#print(el[0])
			if el[0] == fwdrecord.id:

				#open rev index file and seach index seq for matching index name
				with open("tempdir2/index_rev.uniq.renamed.fasta") as indexes_rev:
					for revrecord in SeqIO.parse(indexes_rev, "fasta"):
						#print(el[1])
						if el[1] == revrecord.id:

							#open paired-indexes file
							with open(indexes_file) as indexes_all:
								for indexfilerecord in SeqIO.parse(indexes_all, "fasta"):
									
									index_part = indexfilerecord.seq.split("...")
									# index_part[0] -> forward index seq
									# index_part[1] -> reverse index seq

									if fwdrecord.seq in index_part[0]:
										if revrecord.seq in index_part[1]:
											#print(indexfilerecord.id + ", with indexes " + indexfilerecord.seq)
											sampl_nameR1 = "demultiplex_out/" + str(indexfilerecord.id) + ".R1." + str(extension)
											sampl_nameR2 = "demultiplex_out/" + str(indexfilerecord.id) + ".R2." + str(extension)

											line_R1 = "demultiplex_out/" + line.strip() + ".R1." + str(extension)
											line_R2 = "demultiplex_out/" + line.strip() + ".R2." + str(extension)

											#rename files
											if os.path.isfile(line_R1):
												os.rename(line_R1, sampl_nameR1)
											if os.path.isfile(line_R2):
												os.rename(line_R2, sampl_nameR2)
print("python module finished")
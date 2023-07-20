#!/usr/bin/env Rscript

#DADA2 sequence classifier (module for taxonomy_dada2.sh).
 #input is fasta file

#load dada2
library("dada2")
library("seqinr")

#load env variables
readType = Sys.getenv('readType')
fileFormat = Sys.getenv('fileFormat')
dataFormat = Sys.getenv('dataFormat')
workingDir = Sys.getenv('workingDir')

#check for output dir and delete if needed
if (dir.exists("/input/taxonomy_out.dada2")) {
    unlink("/input/taxonomy_out.dada2", recursive=TRUE)
}
#create output dir
path_results = "/input/taxonomy_out.dada2"
dir.create(path_results)

#load environment variables
database = Sys.getenv('dada2_database')
database = gsub("\\\\", "/", database) #replace backslashes \ in the database path
database = paste("/extraFiles", basename(database), sep = "/")
minBoot = as.integer(Sys.getenv('minBoot'))
tryRC = Sys.getenv('tryRC')
print(database)

#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (tryRC == "false" || tryRC == "FALSE"){
    tryRC = FALSE
}
if (tryRC == "true" || tryRC == "TRUE"){
    tryRC = TRUE
}

#load sequences
if (file.exists("ASVs_lenFilt.fasta") == TRUE && file.exists("ASVs_collapsed.fasta") == TRUE) {
    seqs_file = list.files(file.path(getwd()), pattern = "ASVs_lenFilt.fasta")
} else {
    seqs_file = list.files(file.path(workingDir), pattern = fileFormat)
}

print(seqs_file)


fasta = read.fasta(seqs_file, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE, seqonly = FALSE)
seq_names = getName(fasta)
seqs = unlist(getSequence(fasta, as.string = TRUE))
#Print no of ASVs
paste("Number of sequences = ", length(seq_names))


#assign taxonomy
set.seed(100)
tax = assignTaxonomy(fasta, database, multithread = TRUE, minBoot = minBoot, tryRC = tryRC, outputBootstraps = TRUE)
#add sequence names to tax table
tax2 = cbind(row.names(tax$tax), tax$tax, tax$boot)
colnames(tax2)[1] = "Sequence"
#write taxonomy
write.table(tax2, file.path(path_results, "taxonomy.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

#DONE, proceed with taxonomy_dada2.sh to clean up make readme

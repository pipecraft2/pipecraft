#!/usr/bin/env Rscript

#DADA2 sequence classifier (module for taxonomy_dada2.sh).
 #input is fasta file

#load dada2
library("dada2")
library("seqinr")

#print DADA2 version
cat("DADA2 version = ", base::toString(packageVersion("dada2")), "\n")

#load env variables
workingDir = Sys.getenv('workingDir')
path_results = Sys.getenv('output_dir')

#check for output dir and delete if needed
if (dir.exists(path_results)) {
    unlink(path_results, recursive = TRUE)
}
#create output dir
dir.create(path_results)

#load environment variables
database = Sys.getenv('dada2_database')
database = gsub("\\\\", "/", database) #replace backslashes \ in the database path
database = paste("/extraFiles", basename(database), sep = "/")
fasta_file = Sys.getenv('fasta_file')
fasta_file = gsub("\\\\", "/", fasta_file) #replace backslashes \ in the fasta file path
fasta_file = paste("/extraFiles2", basename(fasta_file), sep = "/")
minBoot = as.integer(Sys.getenv('minBoot'))
tryRC = Sys.getenv('tryRC')

#"FALSE" or "TRUE" to FALSE or TRUE for dada2
if (tryRC == "false" || tryRC == "FALSE"){
    tryRC = FALSE
}
if (tryRC == "true" || tryRC == "TRUE"){
    tryRC = TRUE
}

#log
cat(";; input = ", fasta_file, "\n")
cat(";; database file = ", database, "\n")

#read.fasta
fasta = read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE, seqonly = FALSE)
seq_names = getName(fasta)
seqs = unlist(getSequence(fasta, as.string = TRUE))
#Print no of sequences in the input file
cat(";; Number of sequences = ", length(seq_names))

#assign taxonomy
set.seed(100)
cat(";; Assigning taxonomy with assignTaxonomy function ")
tax = assignTaxonomy(seqs, database, multithread = TRUE, minBoot = minBoot, tryRC = tryRC, outputBootstraps = TRUE)

#add sequence names to tax table
tax2 = cbind(rownames(tax$tax), tax$tax, tax$boot)
rownames(tax2) = seq_names
colnames(tax2)[1] = "Sequence"
#write taxonomy
write.csv(tax2, file.path(path_results, "taxonomy.csv"), row.names = TRUE, quote = FALSE)
cat(";; DONE ")
#DONE, proceed with taxonomy_dada2.sh to clean up make readme

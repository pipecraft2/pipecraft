#!/usr/bin/Rscript

# Generate an OTU table based on the clustered DADA2 ASVs (*.uc file).
# Inputs = ASV_table and clusters.uc file
# Output = OTU_table.txt

library(data.table)

## Required inputs
inp_ASVTAB <- "ASV_table.txt"   # ASV table file (tab delimited). 
                                # First column header = "ASV", next columns = samples.
                                # NO SEQUENCE COLUMN! 
                                # NO SIZE ANNOTATION of ASVs.
inp_UC     <- "OTUs.uc"         # output from `vsearch --cluster_size ... -uc OTU.uc`

#output
out_OTUTAB <- "OTU_table.txt"   # resulting OTU table name 

################################
## Load input data - ASV table
ASVTAB <- fread(file = inp_ASVTAB, header = TRUE, sep = "\t")

#drop 2nd col which has seqs (PipeCraft output table) --- HAS to done in bash
ASVTAB[[2]] <- NULL

## Load input data - UC mapping file
UC <- fread(file = inp_UC, header = FALSE, sep = "\t")
UC <- UC[ V1 != "S" ]
UC[, ASV := tstrsplit(V9, ";", keep = 1) ]
UC[, OTU := tstrsplit(V10, ";", keep = 1) ]
UC[V1 == "C", OTU := ASV ]
UC <- UC[, .(ASV, OTU)]

## Convert ASV table to long format
  #V1 is the col.names for ASVs column (cuz input had not col name)
ASV <- melt(data = ASVTAB,
	id.vars = "V1",
	variable.name = "SampleID", value.name = "Abundance")
#Rename V1 col to "ASV"
colnames(ASV)[1] = "ASV"
ASV <- ASV[ Abundance > 0 ]

## Add OTU IDs
ASV <- merge(x = ASV, y = UC, by = "ASV", all.x = TRUE)

## Summarize
OTU <- ASV[ , .(Abundance = sum(Abundance, na.rm = TRUE)), by = c("SampleID", "OTU")]

## Reshape to wide format
RES <- dcast(data = ASV,
  formula = OTU ~ SampleID,
  value.var = "Abundance",
  fun.aggregate = sum, fill = 0)

## Export OTU table
  # OTU names correspont to most abundant ASV in an OTU
fwrite(x = RES, file = out_OTUTAB, sep = "\t")

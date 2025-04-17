#!/usr/bin/env Rscript

### Parse MetaWorks output "results.csv"

# load workingDir
output_dir = Sys.getenv('output_dir')
cat(";; output_dir = ", output_dir)

# Read the MetaWorks results
results <- read.csv(file.path(output_dir, "results.csv"), header = TRUE)

cat(";; Parsing MetaWorks results; ASV table ")
# Reshape the results into an 'ASV table'
ASV_table <- reshape2::dcast(results, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# replace _ with " " in the results.csv data frame
results_sub <- data.frame(lapply(results, function(x) {gsub("_", " ", x)}))
#convert ESVsize back to integer
results_sub$ESVsize = as.integer(results_sub$ESVsize)
#reshape fun table with taxonomy

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! do the edit: if pseudogene = nom then ORFseq = EMPTY!! 
    # "!"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" ! 
full_table = reshape2::dcast(results_sub, SampleName ~ GlobalESV + ORFseq +
                                Root + rBP + SuperKingdom + skBP + Kingdom +
                                kBP + Phylum + pBP + Class + cBP +
                                Order + oBP + Family + fBP + Genus +
                                gBP + Species + sBP,
                                value.var = "ESVsize", fun.aggregate = sum)

# set ',' as a separator for ASV + Seq + Tax + bootstrap values
names(full_table) = stringr::str_replace_all(names(full_table), "_", ",")
#transpose full table
t_full_table = t(full_table)

cat(";; Parsing MetaWorks results; taxonomy table ")
#take only tax info and discard first row (=SampleName)
tax = rownames(t_full_table)[2:nrow(as.data.frame(t_full_table))] 
#convert to data.frame for tidyr
df.tax = as.data.frame(tax, col.names = "tax")
#split tax columns
tax_table = tidyr::separate(df.tax, 
                            col = "tax", 
                            sep = ",", 
                            into = c("GlobalESV", "ORFseq",
                                    "Root", "rBP",
                                    "SuperKingdom", "skBP",
                                    "Kingdom", "kBP",
                                    "Phylum", "pBP",
                                    "Class", "cBP",
                                    "Order", "oBP",
                                    "Family", "fBP",
                                    "Genus", "gBP",
                                    "Species", "sBP", NA), 
                                    fill = "right")

# Write tables to output and order records
cat(";; Writing output tables")
ASV_table_t = t(ASV_table)
write.table(ASV_table_t[order(rownames(ASV_table_t)),],
"final_ESV_table.txt", row.names = TRUE, col.names = FALSE)

write.table(tax_table[order(tax_table$GlobalESV),],
"final_tax_table.txt", row.names = FALSE, col.names = TRUE)

cat(";; DONE")

#!/usr/bin/env Rscript

#Filter DADA2 ASV table; collapse no mismatch and/or filer by ASV length. For DADA2 full workflow.

#load variables
table_rds = "/input/ASVs_out.dada2/ASVs_table.denoised.nochim.rds"
collapseNoMismatch = Sys.getenv('collapseNoMismatch')
minOverlap = as.numeric(Sys.getenv('minOverlap'))
vec = Sys.getenv('vec')
len_filt = as.numeric(Sys.getenv('by_length'))

#ASV table path
path_out="/input/ASVs_out.dada2/filtered"

#load dada2
library('dada2')
print("### Filtering the ASV table")
print("")

##########################################
#Combine together sequences that are identical up to shifts and/or length.
if (collapseNoMismatch == "true") {
    print("Combining together ASVs that are identical")
    
    # Load ASVs table from the DADA2 workflow
    table_in = readRDS(table_rds)

    ASV_tab_collapsed = collapseNoMismatch(table_in, minOverlap = minOverlap, vec = vec)
    saveRDS(ASV_tab_collapsed, file.path(path_out, "ASV_tab_collapsed.rds"))

    #Print ASV count
    print(paste0("total no. of input ASVs = ", dim(table_in)[2]))
    print(paste0("no. of ASVs in a collapsed table = ", dim(ASV_tab_collapsed)[2]))
    print("")

    ###format and save ASV table and ASVs.fasta
    #sequence headers
    asv_seqs = colnames(ASV_tab_collapsed)
    asv_size = colSums(ASV_tab_collapsed)
    asv_headers = vector(dim(ASV_tab_collapsed)[2], mode="character")
    for (i in 1:dim(ASV_tab_collapsed)[2]) {
        asv_headers[i] = paste(">ASV", i, ";size=", asv_size[i], sep="")
    }
    #transpose sequence table
    tASV_tab_collapsed = t(ASV_tab_collapsed)
    #add sequences to 1st column
    tASV_tab_collapsed = cbind(row.names(tASV_tab_collapsed), tASV_tab_collapsed)
    colnames(tASV_tab_collapsed)[1] = "Sequence"
    #row names as sequence headers
    row.names(tASV_tab_collapsed) = sub(">", "", asv_headers)
    #write ASVs.fasta to path_out
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, file.path(path_out, "ASVs_collapsed.fasta"))
    #write ASVs table to path_out
    write.table(tASV_tab_collapsed, file.path(path_out, "ASVs_table_collapsed.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

    #input for length filtering, if != 0
    table_in = ASV_tab_collapsed

    #remove datasets from the environment
    rm(tASV_tab_collapsed)
    rm(asv_seqs)
    rm(asv_headers)
    rm(asv_fasta)
}


##########################################

##########################################
# Filter ASVs based on length
if (collapseNoMismatch != "true") {
    table_in = readRDS(table_rds)
}

if (len_filt != 0) {
    print(paste0("Filtering ASVs by length; discarding ASVs with sequence length < ", len_filt))
    a = c()
    short_ASV_count = 0
    for (i in 1:dim(table_in)[2]) {
        #if len < len_filt 
        x = as.numeric(nchar(colnames(table_in)[i])) #ASV seq len
        if (x < len_filt) {
            a = append(a, i)
            short_ASV_count = short_ASV_count + 1 
        }
    }

    #Write output file indicationg that no ASVs were filtered out based on this length threshold
    if (is.null(a) == TRUE) {
        write(a, file.path(path_out, "a.txt"))
    }

    #a = NULL if all ASVs were kept. Proceed if there are some ASVs to be removed; i.e. a != NULL
    if (is.null(a) != TRUE) {
        ASV_tab_lenFilt = table_in[,-c(a)] #remove columns, i.e ASVs with short seqs

        #print ASV count
        print(paste0("< ", len_filt, " bp ASVs = ", short_ASV_count))
        print(paste0("no. of ASVs in a length filtered table = ", dim(ASV_tab_lenFilt)[2]))
        print("")

        #Proceed if NOT all ASVs were removed by length filtering
        if (dim(ASV_tab_lenFilt)[2] != 0) { 
            ###format and save ASV table and ASVs.fasta
            #sequence headers with size
            asv_seqs = colnames(ASV_tab_lenFilt)
            asv_size = colSums(ASV_tab_lenFilt)
            asv_headers = vector(dim(ASV_tab_lenFilt)[2], mode="character")
            for (i in 1:dim(ASV_tab_lenFilt)[2]) {
                asv_headers[i] = paste(">ASV", i, ";size=", asv_size[i], sep="")
            }
            #transpose sequence table
            tASV_tab_lenFilt = t(ASV_tab_lenFilt)
            #add sequences to 1st column
            tASV_tab_lenFilt = cbind(row.names(tASV_tab_lenFilt), tASV_tab_lenFilt)
            colnames(tASV_tab_lenFilt)[1] = "Sequence"
            #row names as sequence headers
            row.names(tASV_tab_lenFilt) = sub(">", "", asv_headers)

            #write ASVs.fasta to path_out
            asv_fasta <- c(rbind(asv_headers, asv_seqs))
            write(asv_fasta, file.path(path_out, "ASVs_lenFilt.fasta"))
            #write ASVs table to path_out
            write.table(tASV_tab_lenFilt, file.path(path_out, "ASV_table_lenFilt.txt"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
        } else {
            print(paste0("NO ASVs remained after length filtering; ", len_filt, " bp"))
        }
    }

}
##########################################

#DONE 


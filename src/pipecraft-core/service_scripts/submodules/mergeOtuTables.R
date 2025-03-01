#!/usr/bin/env Rscript

####################################################################
# function to merge multiple OTU tables
# input OTU table format:
# OTU_ID    Sequence    Sample1    Sample2    Sample3
# seq1      TGCGTACGTA    10         20         30
# seq2      TGCGTACGTA    15         25         35
# seq3      TGCGTACGTA    20         30         40
# seq4      TGCGTACGTA    25         35         45
####################################################################
mergeOtuTables <- function(tables = NULL, merge_method = "error") {
    # Input validation
    if (is.null(tables) || length(tables) < 2) {
        stop("At least two OTU tables are required")
    }
    if (!merge_method %in% c("error", "sum")) {
        stop("merge_method must be either 'error' or 'sum'")
    }

    # Read tables if file paths are provided
    if (is.character(tables)) {
        cat(";; Reading", length(tables), "tables...\n")  # Progress message
        tables <- lapply(tables, function(x) {
            df <- read.table(x, header = TRUE, sep = "\t", check.names = FALSE)
            # Check for empty table
            if (nrow(df) == 0) {
                stop("Empty table found in: ", x)
            }
            # Check minimum required columns
            if (ncol(df) < 3) {
                stop("Table must have at least 3 columns (OTU, Sequence, and one sample): ", x)
            }
            # Store OTU and Sequence columns
            meta_cols <- df[, 1:2]
            # Convert abundance columns to numeric
            abundance_cols <- df[, -(1:2), drop = FALSE]
            list(meta = meta_cols, abundance = abundance_cols)
        })
    }

    # Validate table format and convert to numeric
    for (i in seq_along(tables)) {
        if (!is.list(tables[[i]])) {
            stop("Table ", i, " is not properly formatted")
        }
        # Convert all abundance columns to numeric
        tables[[i]]$abundance <- as.data.frame(lapply(tables[[i]]$abundance, function(x) {
            as.numeric(as.character(x))
        }))
        # Check if conversion produced any NAs
        if (any(is.na(tables[[i]]$abundance))) {
            stop("Table ", i, " contains non-numeric values in abundance columns")
        }
    }

    # Check for duplicate sample names if merge_method is "error"
    if (merge_method == "error") {
        all_samples <- lapply(tables, function(x) colnames(x$abundance))
        for (i in 1:(length(all_samples)-1)) {
            for (j in (i+1):length(all_samples)) {
                common_samples <- intersect(all_samples[[i]], all_samples[[j]])
                if (length(common_samples) > 0) {
                    stop("Duplicate sample names found between tables ", i, " and ", j, 
                         ": ", paste(common_samples, collapse=", "))
                }
            }
        }
    }

    # Get all unique OTUs and their sequences
    all_meta <- do.call(rbind, lapply(tables, function(x) x$meta))
    unique_meta <- unique(all_meta)
    
    # Get all unique sample names
    sample_names <- unique(unlist(lapply(tables, function(x) colnames(x$abundance))))

    # Initialize result matrix with numeric type
    merged_abundance <- matrix(0, 
                             nrow = nrow(unique_meta), 
                             ncol = length(sample_names),
                             dimnames = list(NULL, sample_names))

    # Merge abundance values
    for (tab in tables) {
        tab_matrix <- as.matrix(tab$abundance)
        storage.mode(tab_matrix) <- "numeric"
        storage.mode(merged_abundance) <- "numeric"
        
        # Match OTUs between current table and merged result
        matches <- match(paste(tab$meta[,1], tab$meta[,2]), 
                        paste(unique_meta[,1], unique_meta[,2]))
        
        # For overlapping samples and OTUs, sum the values
        merged_abundance[matches, colnames(tab$abundance)] <- 
            merged_abundance[matches, colnames(tab$abundance)] + tab_matrix
    }

    # Convert to final data frame
    merged_table <- data.frame(unique_meta, merged_abundance, check.names = FALSE)
    colnames(merged_table)[1:2] <- c("OTU", "Sequence")

    # Sort by total abundance
    row_sums <- rowSums(merged_table[,-(1:2), drop = FALSE])
    merged_table <- merged_table[order(row_sums, decreasing = TRUE),]

    return(merged_table)
}
####################################################################
merge_method = "sum" # "error" or "sum". "error" will stop if duplicate sample names are found. "sum" will sum samples with same name.

input_tables <- Sys.getenv("output_feature_table")
output_dir <- Sys.getenv("output_dir")

# get input tables (txt files)
cat(";; output_feature_table: ", input_tables, "\n")

# Split the comma-separated string into a vector of file paths
input_table_paths <- strsplit(input_tables, ",")[[1]]
input_table_paths <- trimws(input_table_paths)  # Remove any whitespace

# Check if files exist
missing_files <- input_table_paths[!file.exists(input_table_paths)]
if (length(missing_files) > 0) {
    stop("Error: Following files not found: ", paste(missing_files, collapse=", "))
}

# Merge tables
cat(";; Merging OTU tables...\n")
merged_result <- mergeOtuTables(tables = input_table_paths, merge_method = merge_method)

# Write output
cat(";; Writing merged OTU table...\n")
write.table(merged_result, 
            file = file.path(output_dir, "OTU_table.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat(";; Done!\n")

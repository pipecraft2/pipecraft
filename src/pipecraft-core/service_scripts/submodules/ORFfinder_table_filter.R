#!/usr/bin/env Rscript
### Filter dataset based on ORF-finder results to exclude putative NUMTs

# arguments (paths are provided by ORFfinder.sh)
args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

feature_table_file <- get_arg_value("--table_file")
discard_file <- get_arg_value("--notORFs_list")
outdir <- get_arg_value("--output_dir", default = ".")

if (is.null(feature_table_file) || feature_table_file == "") {
  stop("Missing required argument: --table_file", call. = FALSE)
}
if (is.null(discard_file) || discard_file == "") {
  stop("Missing required argument: --notORFs_list", call. = FALSE)
}

if (!file.exists(feature_table_file)) {
  stop(paste0("Table file not found: ", feature_table_file), call. = FALSE)
}
if (!file.exists(discard_file)) {
  stop(paste0("notORFs list not found: ", discard_file), call. = FALSE)
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Read the list of ASVs that should be discarded based on ORFfinder output
discard = read.table(discard_file, header = FALSE)

# Read the ASV table
feature_table = read.table(feature_table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
n_input = nrow(feature_table)

# Use first column as rownames (ASV IDs), then drop it from the table
rownames(feature_table) <- feature_table[[1]]
feature_table <- feature_table[, -1, drop = FALSE]

# Drop sequence column if present
if ("Sequence" %in% colnames(feature_table)) {
  feature_table <- feature_table[, colnames(feature_table) != "Sequence", drop = FALSE]
}

# Remove non-ORFs from feature table
discard_ids <- unique(trimws(as.character(discard[[1]])))
discard_ids <- discard_ids[!is.na(discard_ids) & discard_ids != ""]
feature_table <- feature_table[!rownames(feature_table) %in% discard_ids, , drop = FALSE]

# Further filter feature_table to remove OTUs/ASVs with zero total reads
if (nrow(feature_table) > 0) {
  row_totals <- rowSums(feature_table, na.rm = TRUE)
  feature_table <- feature_table[row_totals > 0, , drop = FALSE]
}

# Further filter feature_table to remove samples (columns) with zero total reads
if (nrow(feature_table) > 0 && ncol(feature_table) > 0) {
  col_totals <- colSums(feature_table, na.rm = TRUE)
  feature_table <- feature_table[, col_totals > 0, drop = FALSE]
}

### Write outputs
# Output table filename: input basename with "_ORFs" before extension
output_table_filename <- function(f) {
  paste0(tools::file_path_sans_ext(basename(f)), ".ORFs.", tools::file_ext(f))
}

feature_table_out <- data.frame(
  Feature = rownames(feature_table),
  feature_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(feature_table_out, file.path(outdir, output_table_filename(feature_table_file)),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

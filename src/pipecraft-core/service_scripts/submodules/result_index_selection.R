#!/usr/bin/env Rscript

## submodule for metaMATE, when running 'find_and_dump'
# get the results_index from the metamate_out/results.csv file that corresponds to the specified $NA_abund_thresh
# NA_abund_thresh is the allowed abundance threshold of non-validated (putative artefactual) OTUs/ASVs in the filtered dataset.
 
cat("Running result_index_selection.R\n")
## read results.csv
output_dir = Sys.getenv('output_dir')
find_results = read.csv(file.path(output_dir, "results.csv"))
# get the abundance_filt from the environment variable
abundance_filt = Sys.getenv('abundance_filt')

### filter results if abundance_filt is FALSE
if (abundance_filt != "true"){
    result_index = "0" # get first result_index (library_n = 0)
    write(result_index, file.path(output_dir, "selected_result_index.txt"))
    cat("abundance_filt = ", abundance_filt, "\n")
    cat("Selected result_index = 0\n")
}

### filter results based on NA_abund_thresh 
if (abundance_filt == "true"){
    NA_abund_thresh = as.numeric(Sys.getenv('NA_abund_thresh'))
    filtered_data = find_results[find_results$nonauthentic_retained_estimate_p <= NA_abund_thresh, ] 
    cat("abundance_filt = ", abundance_filt, "\n")
    # if no results correspond with the NA_abund_thresh, then get the next best
        # else, just select the result_index that corresponds to NA_abund_thresh with highest accuracy_score
    if (nrow(filtered_data) == 0) {
        cat("\n no results correspond with the NA_abund_thresh of", NA_abund_thresh, "; getting the next best setting\n")
        next_best = min(find_results$nonauthentic_retained_estimate_p)
        filtered_data = find_results[find_results$nonauthentic_retained_estimate_p <= next_best, ] 
        # sort based on accuracy_score
        sorted_filtered = filtered_data[order(-filtered_data$accuracy_score), ]
        # get the result with the highest accuracy_score
        metamate_selected_threshold = sorted_filtered[1,]
        write.csv(metamate_selected_threshold, file.path(output_dir, "next_best_set.csv"), quote = F)
        # the result_index of the NA_abund_thresh with the highest accuracy_score
        result_index = metamate_selected_threshold[,1]
        write(result_index, file.path(output_dir, "selected_result_index.txt"))
    } else {
        # sort based on accuracy_score
        sorted_filtered = filtered_data[order(-filtered_data$accuracy_score), ]
        # get the result with the highest accuracy_score
        metamate_selected_threshold = sorted_filtered[1,]
        # the result_index of the NA_abund_thresh with the highest accuracy_score
        result_index = metamate_selected_threshold[,1]
        write(result_index, file.path(output_dir, "selected_result_index.txt"))
    }
    cat("Selected result_index = ", result_index, "\n")
}
cat("DONE, result_index_selection.R\n")

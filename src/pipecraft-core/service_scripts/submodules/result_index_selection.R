#!/usr/bin/env Rscript

## submodule for metaMATE, when running 'find_and_dump'
# get the results_index from the metamate_out/results.csv file that corresponds to the specified $NA_abund_thresh
# NA_abund_thresh is the allowed abundance threshold of non-validated (putative artefactual) OTUs/ASVs in the filtered dataset.
 
## read results.csv
NA_abund_thresh = Sys.getenv('NA_abund_thresh')
output_dir = Sys.getenv('output_dir')
find_results = read.csv(file.path(output_dir, "results.csv"))

## count ASVs, vaASNs and NUMTs
#ASVs_total = find_results$asvs_total[1]
#vaASVs_total = find_results$verifiedauthentic_total_observed[1]
#NUMTs_total = find_results$verifiednonauthentic_total_observed[1]

## filter results based on NA_abund_thresh 
filtered_data = find_results[find_results$nonauthentic_retained_estimate_p <= NA_abund_thresh, ] 
# sort based on accuracy_score
sorted_filtered = filtered_data[order(-filtered_data$accuracy_score), ]
# get the result with the highest accuracy_score
metamate_selected_threshold <- sorted_filtered[1,]
# the result_index of the NA_abund_thresh with the highest accuracy_score
result_index = metamate_selected_threshold[,1]
write(result_index, file.path(output_dir, "selected_result_index.txt"))

# if no results correspond with the NA_abund_thresh, then get the next best !!!

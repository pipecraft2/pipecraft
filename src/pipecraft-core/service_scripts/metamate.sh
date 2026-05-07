#!/bin/bash

# metaMATE v0.5.4 (https://github.com/tjcreedy/metamate) software for 
 #  removel of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data

#load variables
filter_mode=${filter_mode:-$per-sample}   # global, per-sample
specifications=${specifications}   # file
reference_seqs=${reference_seqs}   # file
reference_seqs2=${reference_seqs2}   # file (optional)
table=${table}                     # file 
rep_seqs=${rep_seqs}               # file
genetic_code=${genetic_code}       # integer
length=${length}                   # integer
result_index=${result_index}       # integer --_> DELETE!
abundance_filt=${abundance_filt}   # boolean
NA_abund_thresh=${NA_abund_thresh} # float
bases_variation=${bases_variation} # integer
taxgroups=${taxgroups}             # file (optional)
cores=${cores}                     # integer
percentile=${percentile}           # float (0-1), filter-adaptive
criteria=${criteria}               # verified_removed | estimated_removed, filter-adaptive
realign=${realign}                 # boolean (metaMATE --realign)
keeptemporaryfiles=${keeptemporaryfiles} # boolean (metaMATE --keeptemporaryfiles)
otu_mode=${otu_mode}               # boolean
uc=${uc}                           # file (.uc)
otu_fasta=${otu_fasta}             # file (fasta)
otu_table=${otu_table}             # file (csv/tsv)

# Function to convert a UI-provided host file path into the container-visible path.
# Strategy: prefer `/input/<basename>` (working directory mount), otherwise search `/extraFiles*`
container_file() {
    local p="$1"
    if [[ "$p" == "undefined" ]] || [[ -z "$p" ]]; then
        echo ""
        return
    fi

    local base
    base="$(basename "$p")"

    # Prefer the selected working directory mount.
    if [[ -f "/input/$base" ]]; then
        echo "/input/$base"
        return
    fi

    # Fall back to any extraFiles mount (slot numbers depend on input order/visibility).
    local d
    for d in /extraFiles*; do
        if [[ -d "$d" ]] && [[ -f "$d/$base" ]]; then
            echo "$d/$base"
            return
        fi
    done

    # Default to /input for clearer error messages downstream.
    echo "/input/$base"
}

source /scripts/submodules/framework.functions.sh

# Specifications file handling
# Rule:
# - If the user does NOT provide a specifications file, always use the bundled
#   image default at `/default_specs/specifications0.txt` (created in Dockerfile).
# - If the user DOES provide a file, use it (mapped into the container).
if [[ "$abundance_filt" != "true" ]]; then
    specifications="/default_specs/specifications0.txt"
    printf '%s\n' "NOTICE]: abundance_filt is not 'true' — forcing default metaMATE specifications: $specifications" >&2
else
    if [[ -z "$specifications" ]] || [[ "$specifications" == "undefined" ]]; then
        specifications="/default_specs/specifications.txt"
        printf '%s\n' "NOTICE]: no specifications file provided — using bundled default: $specifications" >&2
    else
        specifications=$(container_file "$specifications")
        printf '%s\n' "NOTICE]: using user-provided metaMATE specifications: $specifications" >&2
    fi
fi

# Reference seqs (database) handling
if [[ "${reference_seqs}" != "undefined" ]] && [[ -n "${reference_seqs}" ]]; then
    reference_seqs=$(container_file "$reference_seqs")
fi

# Reference seqs2 handling
if [[ $reference_seqs2 != "undefined" ]]; then
    reference_seqs2=$(container_file "$reference_seqs2")
fi

# Feature table handling
if [[ "${table}" != "undefined" ]] && [[ -n "${table}" ]]; then
    table=$(container_file "$table")
fi
# check "Sequence" column in the table file
if [[ -n "$table" ]] && [[ -f "$table" ]] && grep -q "Sequence" "$table"; then
    printf '%s\n' "WARNIG]: table file contains a 'Sequence' column. Removing this column before running metaMATE." >&2
    awk 'BEGIN {FS=OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i=="Sequence") col=i} \
        {for (i=1; i<=NF; i++) if (i!=col) printf "%s%s", $i, (i==NF?ORS:OFS)}' $table > $table.temp
    table=$table.temp
fi

# Rep seqs handling
rep_seqs=$(container_file "$rep_seqs")

if [[ -z "$rep_seqs" ]] || [[ ! -f "$rep_seqs" ]]; then
    printf '%s\n' "ERROR]: rep_seqs is missing or not found in container mounts. Put it in the PipeCraft working directory (/input) or ensure it is selectable as a file input (mounted under /extraFiles*)." >&2
    end_process
fi
if [[ -z "$table" ]] || [[ ! -f "$table" ]]; then
    printf '%s\n' "ERROR]: ASV/OTU table is missing or not found in container mounts. Put it in the working directory (/input) or ensure it is selectable as a file input (mounted under /extraFiles*)." >&2
    end_process
fi

# Taxgroups handling
if [[ $taxgroups != "undefined" ]]; then
    taxgroups=$(container_file "$taxgroups")
    taxgroups=$"--taxgroups $taxgroups"
else 
    taxgroups=$""
fi

# OTU mode handling
otu_args=$""
if [[ $otu_mode == "true" ]]; then
    if [[ $uc == "undefined" ]] || [[ $otu_fasta == "undefined" ]] || [[ $otu_table == "undefined" ]]; then
        printf '%s\n' "ERROR]: OTU mode is enabled, but one or more required files are missing (uc / otu_fasta / otu_table).
>Quitting" >&2
        end_process
    fi

    uc=$(container_file "$uc")
    otu_fasta=$(container_file "$otu_fasta")
    otu_table=$(container_file "$otu_table")

    otu_args=$"--uc $uc --otu_fasta $otu_fasta --otu_table $otu_table"
fi

# Optional metaMATE flags (only apply to find/filter-adaptive runs)
metamate_extra_args=$""
if [[ "$realign" == "true" ]]; then
    metamate_extra_args=$"$metamate_extra_args --realign"
fi
if [[ "$keeptemporaryfiles" == "true" ]]; then
    metamate_extra_args=$"$metamate_extra_args --keeptemporaryfiles"
fi

printf "\n specifications file = $specifications\n"
printf "reference seqs file = $reference_seqs\n"
if [[ $reference_seqs2 != "undefined" ]]; then
    printf "reference seqs file2 = $reference_seqs2\n"
fi
printf "filter mode = $filter_mode\n"
printf "otu mode = $otu_mode\n"
printf "uc file = $uc\n"
printf "otu fasta file = $otu_fasta\n"
printf "otu table file = $otu_table\n"
printf "table file = $table\n"
printf "rep seqs file = $rep_seqs\n"
printf "cores = $cores\n"
printf "length = $length\n"
printf "result_index = $result_index\n"
printf "abundance_filt = $abundance_filt\n"
printf "NA_abund_thresh = $NA_abund_thresh\n"
printf "bases_variation = $bases_variation\n"
printf "taxgroups = $taxgroups\n"
###############################
#output dir
output_dir=$"/input/metamate_out"
echo "output_dir = $output_dir"


#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)

### quick check of the rep_seqs file
if ! grep -q "^>" $rep_seqs; then
    printf '%s\n' "ERROR]: rep_seqs file is not a FASTA file.
    Please provide a fasta.
    >Quitting" >&2
    end_process
fi

### if perfoming clade binning, then WARNING when processing more than 65,536 ASVs
ASVcount=$(grep -c "^>" $rep_seqs)
if (( $ASVcount > 65536 )); then
    printf '%s\n' "WARNING]: clade binning NOT performed, because the input ASVs limit is 65,536 for that.
    Current input has $ASVcount ASVs."
fi

### quick check that the number of ASVs/OTUs match between $rep_seqs and $table files
num_lines=$(wc -l < $table)
ASVcount_in_table=$((num_lines - 1))
if [[ $ASVcount != $ASVcount_in_table ]]; then
    printf '%s\n' "ERROR]: number of ASVs/OTUs do not match between the specified rep_seqs and table files.
Not corresponding files submitted or the table format is wrong. Expected input table format = samples in COLUMNs and ASVs/OTUs in ROWs.
>Quitting" >&2
    end_process
fi

### check the length distribution in the rep_seqs file
    # error if specified length +- base_variation is not within the length distribution
length_min=$(seqkit stats -T $rep_seqs | awk '{print $6}' | sed -n '2p')
echo "length_min = $length_min"
length_max=$(seqkit stats -T $rep_seqs | awk '{print $8}' | sed -n '2p')
echo "length_max = $length_max"
# calculate acceptable length range based on specified length and variation
range_min=$((length - bases_variation))
range_max=$((length + bases_variation))
# check for overlap
if [[ $length_min -le $range_max && $range_min -le $length_max ]]; then
  echo "Sequence length ranges overlap, ok!"
else
  printf '%s\n' "ERROR]: the sequence lenght ranges in the rep_seqs file are outside the acceptable length range.
    Specified length = $length bp +/- $bases_variation bp (acceptable range: $range_min - $range_max bp)
    Found sequences ranging from $length_min - $length_max bp
    >Quitting" >&2
    end_process
fi

### shuffle sequences in rep_seqs file 
    # for metaMATE alignment software, because if first 1000 seqs are same length, then 
    # it will think the seqs are already aligned and gives error
rep_seqs_shuffled=$(basename $rep_seqs)
checkerror=$(seqkit shuffle --quiet -w 0 $rep_seqs \
                > /input/${rep_seqs_shuffled%.*}_shuffled.fasta 2>&1)
check_app_error

### metaMATE-find
if [[ $filter_mode == "global" ]]; then

    echo "filter_mode = $filter_mode"
    # quick check of the specifications file, has to contain "library" | "total" | "clade" | "taxon"
    if ! grep -q -e "library" -e "total" -e "clade" -e "taxon" $specifications; then
        printf '%s\n' "ERROR]: specifications file seems to be wrong. Does not contain any of the terms (library, total, clade, taxon).
        >Quitting" >&2
        end_process
    fi

    # quick check of the reference_seqs file
    if ! grep -q "^>" $reference_seqs; then
        printf '%s\n' "ERROR]: reference_seqs file is not a FASTA file.
        Please provide a fasta.
        >Quitting" >&2
        end_process
    fi
    if [[ $reference_seqs2 != "undefined" ]]; then
        if ! grep -q "^>" $reference_seqs2; then
            printf '%s\n' "ERROR]: reference_seqs2 file is not a FASTA file.
            Please provide a fasta.
            >Quitting" >&2
            end_process
        fi
    fi

    # remove old $output_dir if exists
    if [[ -d "$output_dir" ]]; then
        rm -rf "$output_dir"
    fi
    mkdir -p "$output_dir"
    # merge reference seqs files; if reference_seqs2 is provided
    if [[ $reference_seqs2 != "undefined" ]]; then
        cat "$reference_seqs" "$reference_seqs2" > "$output_dir/reference_seqs_merged.fasta"
        db1=$(basename $reference_seqs)
        db2=$(basename $reference_seqs2)
        reference_seqs="$output_dir/reference_seqs_merged.fasta"
    else 
        db1=$(basename $reference_seqs)
        db2=$""
    fi

    printf "# Running metaMATE-find\n"
    checkerror=$(metamate find \
        --asvs /input/${rep_seqs_shuffled%.*}_shuffled.fasta \
        --readmap $table \
        --specification $specifications \
        --references $reference_seqs \
        --expectedlength $length \
        --basesvariation $bases_variation \
        --table $genetic_code \
        --threads $cores \
        --output "$output_dir" \
        --overwrite $taxgroups $otu_args $metamate_extra_args 2>&1)
    check_app_error

    # remove temp table file and assign table to its original name
    if [[ -e $table.temp ]]; then
        original_table=${table%.temp}  # Only remove .temp extension
        rm $table.temp
        table=$original_table
    fi

    # remove shuffled fasta file
    if [[ -e /input/${rep_seqs_shuffled%.*}_shuffled.fasta ]]; then
        rm /input/${rep_seqs_shuffled%.*}_shuffled.fasta
    fi
fi

### metaMATE-filter-adaptive (per-sample filtering)
if [[ $filter_mode == "per-sample" ]]; then

    echo "filter_mode = $filter_mode"

    # quick check of the reference_seqs file
    if ! grep -q "^>" $reference_seqs; then
        printf '%s\n' "ERROR]: reference_seqs file is not a FASTA file.
        Please provide a fasta.
        >Quitting" >&2
        end_process
    fi

    # remove old $output_dir if exists
    if [[ -d "$output_dir" ]]; then
        rm -rf "$output_dir"
    fi
    mkdir -p "$output_dir"

    # merge reference seqs files; if reference_seqs2 is provided
    if [[ $reference_seqs2 != "undefined" ]] && [[ -n "$reference_seqs2" ]]; then
        cat "$reference_seqs" "$reference_seqs2" > "$output_dir/reference_seqs_merged.fasta"
        reference_seqs="$output_dir/reference_seqs_merged.fasta"
    fi

    printf "# Running metaMATE-filter-adaptive (per-sample filtering)\n"
    checkerror=$(metamate filter-adaptive \
        --asvs /input/${rep_seqs_shuffled%.*}_shuffled.fasta \
        --readmap $table \
        --references $reference_seqs \
        --expectedlength $length \
        --basesvariation $bases_variation \
        --table $genetic_code \
        --threads $cores \
        --output "$output_dir" \
        --overwrite $taxgroups $otu_args $metamate_extra_args \
        --percentile $percentile \
        --criteria $criteria 2>&1)
    check_app_error

    # remove shuffled fasta file
    if [[ -e /input/${rep_seqs_shuffled%.*}_shuffled.fasta ]]; then
        rm /input/${rep_seqs_shuffled%.*}_shuffled.fasta
    fi
fi

### metaMATE-dump 	
if [[ $filter_mode == "global" ]]; then

    echo "filter_mode = $filter_mode"

    # dump output name
    dump_seqs=$(basename $rep_seqs)

    # metaMATE writes outputs as *resultcache and *results.csv 
    resultcache_path=$(ls -1 "$output_dir"/*resultcache 2>/dev/null | head -n 1)
    results_csv_path=$(ls -1 "$output_dir"/*results.csv 2>/dev/null | head -n 1)

    # check for the presence of required metaMATE-find outputs
    if [[ -d "$output_dir" ]] && [[ -n "$resultcache_path" ]] && [[ -f "$resultcache_path" ]] && \
    [[ -n "$results_csv_path" ]] && [[ -f "$results_csv_path" ]]; then

        printf "# Running metaMATE-dump\n"

        # get the result_index for the specified NA_abund_thresh
        export NA_abund_thresh
        export output_dir
        Rlog=$(Rscript /scripts/submodules/result_index_selection.R 2>&1) 
        echo "$Rlog" > "$output_dir/result_index_selection.log"
        wait

        # read result_index
        read -r result_index < "$output_dir/selected_result_index.txt"
        printf " - selcted result_index = $result_index\n"

        # if no results correspond with the NA_abund_thresh, then get the next best 

        # run metaMATE-dump
        checkerror=$(metamate dump \
        --asvs $rep_seqs \
        --resultcache $resultcache_path \
        --output $output_dir/${dump_seqs%.*}_metaMATE.filt \
        --overwrite \
        --resultindex $result_index $otu_args 2>&1)
        check_app_error

        # generate a list of ASV IDs 
        checkerror=$(seqkit seq -n "$output_dir/${dump_seqs%.*}_metaMATE.filt.fasta" > "$output_dir/${dump_seqs%.*}_metaMATE.filt.list" 2>&1)
        check_app_error

        # filter the ASV table; include only the ASVs that are in ${dump_seqs%.*}_metaMATE.filt.list
        out_table=$(basename $table)
        awk -v var="$output_dir/${dump_seqs%.*}" 'NR==1; NR>1 {print $0 | "grep -Fwf "var"_metaMATE.filt.list"}' "$table" > "$output_dir/${out_table%.*}_metaMATE.filt.txt"

        # If metaMATE wrote any filtered table(s), keep them alongside Pipecraft outputs
        if compgen -G "$output_dir/*table*" > /dev/null; then
            mkdir -p "$output_dir/metamate_tables" || true
            cp -f "$output_dir"/*table* "$output_dir/metamate_tables/" 2>/dev/null || true
        fi
    else 
        printf '%s\n' "ERROR]: cannot find the $output_dir (metaMATE-find output) to start metaMATE-dump.
Please check that the correct working directory is specified or run metaMATE-find first. Do not use metamate_out as working directory for metaMATE-dump.
>Quitting" >&2
    end_process
    fi
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    if [[ -e "$output_dir/result_index_selection.log" ]]; then
        rm "$output_dir/result_index_selection.log"
    fi
fi
#basename of rep_seqs
rep_seqs=$(basename $rep_seqs)

# rm metaMATE formatted ASV table - not needed
if [[ -e "$output_dir/${rep_seqs%.*}_ASVcounts.csv" ]]; then
    rm "$output_dir/${rep_seqs%.*}_ASVcounts.csv"
fi

#Make README.txt file
if (( $ASVcount > 65536 )); then
    warn=$(echo "WARNING]: clade binning NOT performed, because the input ASVs limit is 65,536 for that. Current input has $ASVcount ASVs.")
fi
end=$(date +%s)
runtime=$((end-start))

if [[ $filter_mode == "global" ]]; then
    printf "### Used metaMATE-find to detect putative NUMT and other erroneous sequences.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input parameters:
---------------
- filter_mode: ${filter_mode}
- specifications: $(basename $specifications)
- reference_seqs: $db1 $db2
- table: $(basename ${table%.temp})
- rep_seqs: $(basename $rep_seqs)
- genetic_code: ${genetic_code}
- length: ${length}
- abundance_filt: ${abundance_filt}
- NA_abund_thresh: ${NA_abund_thresh}
- bases_variation: ${bases_variation}
---------------

Files in 'metamate_out' directory:
----------------------------------
# *_aligned.fasta = aligned $rep_seqs
# *_UPGMA.tre     = newick-format UPGMA tree file (if clade binning was specified)
# results.csv     = metaMATE find results file, which synthesises all of the results from 
                    applying all combinations of the specified terms and thresholds to the input ASVs, given the control groups of authentic- and non-authentic-ASVs.
                    see this metaMATE documentation link for more info about the results file: https://github.com/tjcreedy/metamate?tab=readme-ov-file#results-find-only 
# resultcache     = needed for dump. This is a compressed text file containing information on 
                    the ASVs rejected or retained for each of the supplied specification terms and threshold sets of a find run.
# *_control.txt   = a two-column tab-separated table recording all ASVs determined to be 
                    validated-authentic or validated-non-authentic. 
*_clades.csv      = a two-column comma-separated table recording the clade grouping for 
                    each input ASV.

# -> more info about the 'find' outputs: https://github.com/tjcreedy/metamate?tab=readme-ov-file#outputs 

# Check out analyse_results_draft.R provided by the metaMATE developers for generating plots from a metaMATE find. 
  https://github.com/tjcreedy/metamate/blob/main/analyse_results_draft.R

$warn

#################################################
###Third-party applications for this process:
# metaMATE v0.5.4
    #citation: Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A.J., Vogler, A.P. and Emerson, B.C. (2021), Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data. Mol Ecol Resour, 21: 1772-1787. https://doi.org/10.1111/1755-0998.13337
    #https://github.com/tjcreedy/metamate
#################################################" > "$output_dir/README.metaMATE-find.txt"
fi

if [[ $filter_mode == "global" ]]; then
    # count input ASVs
    inASV_count=$(grep -c "^>" $rep_seqs)
    rep_seqs=$(basename $rep_seqs)
    # count output ASVs
    outASV_count=$(grep -c "^>" "$output_dir/${dump_seqs%.*}_metaMATE.filt.fasta")
 
    # count total sequences, skipping header and handling potential Sequence column
    nSeqs=$(awk 'BEGIN{FS=OFS="\t"}
        NR==1 {
            # Find Sequence column if it exists
            for(i=1;i<=NF;i++) {
                if($i=="Sequence") seq_col=i
                else if(i>1) {
                    h=$i
                    gsub(/\r/,"",h)
                    gsub(/\n/,"",h)
                    gsub(/^[[:space:]]+|[[:space:]]+$/,"",h)
                    header[i]=h  # Store sample names (sanitized)
                }
            }
            next
        }
        {
            # Process each row
            for(i=2;i<=NF;i++) {
                if(i!=seq_col) {
                    sample_sums[i]+=$i  # Sum per sample
                    total+=$i           # Overall total
                }
            }
        }
        END{
            # Print per-sample sums
            printf "Sequences per sample:\n"
            for(i=2;i<=NF;i++) {
                if(i!=seq_col) printf "%s: %s\n", header[i], sample_sums[i]
            }
            printf "\nTotal sequences: %s\n", total
        }' "$output_dir/${out_table%.*}_metaMATE.filt.txt")

    echo "$nSeqs" > "$output_dir/sequence_counts.txt"

    end=$(date +%s)
    runtime=$((end-start))
    printf "### Used metaMATE-dump to discard putative NUMT and other erroneous sequences based on the specified threshold from metaMATE-find.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input parameters:
---------------
- filter_mode: ${filter_mode}
- specifications: $(basename $specifications)
- reference_seqs: $(basename $reference_seqs)
- table: $(basename ${table%.temp})
- rep_seqs: $(basename $rep_seqs)
- genetic_code: ${genetic_code}
- length: ${length}
- abundance_filt: ${abundance_filt}
- NA_abund_thresh: ${NA_abund_thresh}
- bases_variation: ${bases_variation}
---------------

-Input ($rep_seqs) contained $inASV_count sequences
-metaMATE filtered output containes $outASV_count sequences

Added files to 'metamate_out' directory:
----------------------------------------
# ${dump_seqs%.*}_metaMATE.filt.fasta = output of metaMATE-dump function. 
                                        Containes $outASV_count sequences.
# ${dump_seqs%.*}_metaMATE.filt.list  = a list of sequence IDs from ${dump_seqs%.*}_metaMATE.filt.fasta
# ${out_table%.*}_metaMATE.filt.txt = an ASV/OTU table containing the ASVs/OTUs 
                                      that are in ${dump_seqs%.*}_metaMATE.filt.fasta file.
# selected_result_index.txt = if present, then this file contains the 
                              selected resultindex for results.csv file for metaMATE-dump

# -> more info about the outputs: https://github.com/tjcreedy/metamate?tab=readme-ov-file#outputs 

#################################################
###Third-party applications for this process:
# metaMATE v0.5.4
    #citation: Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A.J., Vogler, A.P. and Emerson, B.C. (2021), Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data. Mol Ecol Resour, 21: 1772-1787. https://doi.org/10.1111/1755-0998.13337
    #https://github.com/tjcreedy/metamate
#################################################" > "$output_dir/README.metaMATE-dump.txt"
fi

if [[ -e "$output_dir/next_best_set.csv" ]]; then
    sed -i "7i\# next_best_set.csv = contains the next best filtering settings as the metaMATE-find results.csv file did not contain NA_abund_thresh of <= $NA_abund_thresh ('nonauthentic_retained_estimate_p')." "$output_dir/README.metaMATE-dump.txt"
fi


### README for metaMATE filter-adaptive (per-sample) mode
if [[ $filter_mode == "per-sample" ]]; then

    # rename native metaMATE outputs to more informative names
    mv "$output_dir/filter-adaptive_table_filtered.txt" "$output_dir/feature_table.txt"
    mv "$output_dir/filter-adaptive.fasta" "$output_dir/features.fasta"
    # rename "*_control.txt" to passes_and_fails.txt
    control_candidates=("$output_dir"/*_control.txt)
    if [[ ${#control_candidates[@]} -eq 1 ]] && [[ -f "${control_candidates[0]}" ]]; then
        mv "${control_candidates[0]}" "$output_dir/passes_and_fails.txt"
    fi
    # remove the csv files, not needed
    rm -f "$output_dir/filter-adaptive_table.csv"
    rm -f "$output_dir"/*_ASVcounts.csv

    # count input ASVs
    inASV_count=$(grep -c "^>" "$rep_seqs")
    rep_seqs=$(basename "$rep_seqs")

    # count input OTUs
    if [[ $otu_mode == "true" ]]; then
        inOTU_count=$(grep -c "^>" "$otu_fasta")
        otu_fasta=$(basename "$otu_fasta")
        otu_table=$(basename "$otu_table")
        uc=$(basename "$uc")
    fi

    # count output features
    outASV_count=$(grep -c "^>" "$output_dir/features.fasta")
 
    # count total sequences, skipping header and handling potential Sequence column
    nSeqs=$(awk 'BEGIN{FS=OFS="\t"}
        NR==1 {
            # Find Sequence column if it exists
            for(i=1;i<=NF;i++) {
                if($i=="Sequence") seq_col=i
                else if(i>1) {
                    h=$i
                    gsub(/\r/,"",h)
                    gsub(/\n/,"",h)
                    gsub(/^[[:space:]]+|[[:space:]]+$/,"",h)
                    header[i]=h  # Store sample names (sanitized)
                }
            }
            next
        }
        {
            # Process each row
            for(i=2;i<=NF;i++) {
                if(i!=seq_col) {
                    sample_sums[i]+=$i  # Sum per sample
                    total+=$i           # Overall total
                }
            }
        }
        END{
            # Print per-sample sums
            printf "Sequences per sample:\n"
            for(i=2;i<=NF;i++) {
                if(i!=seq_col) printf "%s: %s\n", header[i], sample_sums[i]
            }
            printf "\nTotal sequences: %s\n", total
        }' "$output_dir/feature_table.txt")

    echo "$nSeqs" > "$output_dir/sequence_counts.txt"

    end=$(date +%s)
    runtime=$((end-start))
    
    cat > "$output_dir/README.metaMATE.per-sample.txt" <<EOF
### Used metaMATE filter-adaptive mode to discard putative NUMT and other erroneous sequences.

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Input parameters:
---------------
- filter_mode: ${filter_mode}
- reference_seqs: $(basename $reference_seqs)
- table: $(basename ${table%.temp})
- rep_seqs: $(basename $rep_seqs)
- genetic_code: ${genetic_code}
- length: ${length}
- bases_variation: ${bases_variation}
- percentile: ${percentile}
- criteria: ${criteria}
- otu mode: ${otu_mode}
- uc file: ${uc}
- otu fasta file: ${otu_fasta}
- otu table file: ${otu_table}
---------------
EOF

    if [[ $otu_mode == "true" ]]; then
        printf "\n-Input ($otu_fasta) contained $inOTU_count OTUs.
-metaMATE filtered output containes $outASV_count OTUs.\n" >> "$output_dir/README.metaMATE.per-sample.txt"
    else
        printf "\n-Input ($rep_seqs) contained $inASV_count features.
-metaMATE filtered output containes $outASV_count features.\n" >> "$output_dir/README.metaMATE.per-sample.txt"
    fi

    cat >> "$output_dir/README.metaMATE.per-sample.txt" <<EOF

Added files to 'metamate_out' directory:
----------------------------------------
# features.fasta = metaMATE filtered features. Containes $outASV_count features.
# feature_table.txt = metaMATE filtered feature table.
# passes_and_fails.txt = list of features with refpass, lenghtfails and stopfails.
                            refpass = feature that perfectly matched the reference sequence
                            lenghtfails = feature did not pass the length filter
                            stopfails = feature did not pass the stop codon filter
# filter-adaptive_summary.csv = summary of applied thresholds, authentic and non-authentic features for each sample
# otu_summary.csv [if otu mode = TRUE ] = summary of ASV and OTU_Status (Authentic, Non-Authentic, Unclassified.
                                          Authentic: feature that perfectly matched the reference sequence
                                          Non-Authentic: feature that did not pass the genetic code translation or length filter
                                          Unclassified: feature that was not classified as Authentic or Non-Authentic

# -> more info about the outputs: https://github.com/tjcreedy/metamate?tab=readme-ov-file#outputs

#################################################
###Third-party applications for this process:
# metaMATE v0.5.4
    #citation: Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A.J., Vogler, A.P. and Emerson, B.C. (2021), Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data. Mol Ecol Resour, 21: 1772-1787. https://doi.org/10.1111/1755-0998.13337
    #https://github.com/tjcreedy/metamate
#################################################
EOF
fi

# Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=fasta"
echo "readType=single_end"

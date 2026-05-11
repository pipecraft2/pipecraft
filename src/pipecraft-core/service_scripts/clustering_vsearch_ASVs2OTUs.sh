#!/bin/bash

# ASVs to OTUs with vsearch
#Input = ASV fasta and ASV table
    # ASVs.fasta must be without size annotations; but size annotation are important, current workflow gets those from the ASV table.
    # ASV table format = ASVs in rows, samples in columns; 2nd column must be SEQUENCE column (default PipeCraft2 ASVs workflow output)
#Output = FASTA formated representative OTU sequences and OTU table

# 11.05.2026: support the case where ASV_fasta is a subset of ASV_table.

#############################
###Third-party applications:
# vsearch
#############################
# checking tool versions
vsearch_version=$(vsearch --version 2>&1 | head -n 1 | awk '{print $2}' | sed -e "s/,//g")
printf "# vsearch (version $vsearch_version)\n"

# Source for functions
source /scripts/submodules/framework.functions.sh

#load ASV_fasta
ASV_fasta=${ASV_fasta}
if [[ $ASV_fasta == "undefined" ]]; then
    printf '%s\n' "ERROR]: ASV_fasta not selected. >Quitting" >&2
    end_process
else
    regex='[^/]*$'
    ASV_fasta_temp=$(echo $ASV_fasta | grep -oP "$regex")
    ASV_fasta=$(printf "/extraFiles/$ASV_fasta_temp")
    printf "\n ASV_fasta = $ASV_fasta \n"
fi
extension=$(basename $ASV_fasta | awk 'BEGIN{FS="."}{print $NF}')

#load ASV_table
ASV_table=${ASV_table}
if [[ $ASV_table == "undefined" ]]; then
    printf '%s\n' "ERROR]: ASV_table not selected. >Quitting" >&2
    end_process
else
    regex='[^/]*$'
    ASV_table_temp=$(echo $ASV_table | grep -oP "$regex")
    ASV_table=$(printf "/extraFiles2/$ASV_table_temp")
    printf "\n ASV_table = $ASV_table \n"
fi

#settings
id=$"--id ${similarity_threshold}"          # positive float (0-1)
otutype=$"--${OTU_type}"                    # list: --centroids, --consout
strands=$"--strand ${strands}"              # list: --strand both, --strand plus
remove_singletons=$"${remove_singletons}"   # true/false
seqsort=$"${sequence_sorting}"              # list: --cluster_size or --cluster_fast, --cluster_smallmem
simtype=$"--iddef ${similarity_type}"       # list: --iddef 0; --iddef 1; --iddef 2; --iddef 3; --iddef 4
centroid=$centroid_type                     # list: similarity, abundance
maxaccepts=$"--maxaccepts ${maxaccepts}"    # pos integer
mask=$"--qmask ${mask}"                     # list: --qmask dust, --qmask none
cores=$"--threads ${cores}"                 # pos integer
###############################

#output dir
output_dir=$"/input/ASVs2OTUs_out"
export output_dir

#additional options, if selection != undefined/false
if [[ $seqsort == "size" ]]; then
    seqsort=$"--cluster_size"
elif [[ $seqsort == "length" ]]; then
    seqsort=$"--cluster_fast"
elif [[ $seqsort == "none" ]]; then
    seqsort=$"--cluster_smallmem --usersort"
fi 
if [[ $centroid == "similarity" ]]; then
    centroid_in=$"" 
else
    centroid_in=$"--sizeorder"
fi
if [[ $remove_singletons == "true"  ]]; then
    remove_singletons=$"TRUE"
fi
if [[ $remove_singletons == "false"  ]]; then
    remove_singletons=$"FALSE"
fi

#############################
### Start of the workflow ###
#############################
start_time=$(date)
start=$(date +%s)

echo "output_dir = $output_dir"
if [[ -d $output_dir ]]; then
    rm -rf $output_dir
fi
mkdir $output_dir
if [[ -d "tempdir" ]]; then
    rm -rf tempdir
fi
mkdir -p tempdir

### Check the ASV table
# 2nd column must be the sequence column.
first_seq=$(awk 'NR==2{print $2; exit}' "$ASV_table" | tr -d '\r')
if [[ -z "$first_seq" ]]; then
    printf '%s\n' "ERROR]: Could not read the 2nd column from $ASV_table.
    >Quitting" >&2
    end_process
fi

# Allow standard IUPAC DNA letters (plus '-' for gaps) and require a reasonable length.
if ! printf '%s' "$first_seq" | grep -Eq '^[ACGTRYSWKMBDHVNacgtryswkmbdhvn-]+$' || [[ ${#first_seq} -lt 20 ]]; then
    printf '%s\n' "ERROR]: 2nd column does not look like a DNA sequence in the provided feature table.
    Detected value in row 2, col 2: '$first_seq'
    If your table does not yet contain a sequence column,
    run QuickTools --> Utilities --> add sequences to table.
    >>>Quitting" >&2
    end_process
fi
# Remove size annotations from the ASV table automatically, if present.
# Subsequent steps must use the cleaned table (if created).
ASV_table_orig="$ASV_table"
ASV_table_nosize=""
ASV_table_clean="$ASV_table"
if head "$ASV_table" | grep -q ";size="; then
    printf "\nFound ';size=' annotations in %s. Removing them automatically.\n" "$ASV_table"
    ASV_table_nosize="tempdir/$(basename "$ASV_table" | sed -E 's/\.[^.]+$//').nosize.txt"
    sed -E 's/;size=[0-9]+//g' "$ASV_table" > "$ASV_table_nosize"
    ASV_table_clean="$ASV_table_nosize"
fi
ASV_table="$ASV_table_clean"
export ASV_table

### Get size annotations for ASVs for clustering (making ASVs.size.fasta from the ASV table where 2nd column is SEQUENCE).
    # Sum each ROW in the table and print ">" + 1st col + row sum as ";size=" + 2nd col (sequence)
    # NR>1 skips the first row, which are the col names.
fastasize=$(basename "$ASV_fasta" | awk 'BEGIN{FS=OFS="."}NF{NF -=1}1')
awk 'NR>1{for(i=3;i<=NF;i++) t+=$i; print ">"$1";size="t"\n"$2; t=0}' "$ASV_table" > "$output_dir/$fastasize.size.fasta"

# Check if $fastasize.size.fasta and $ASV_fasta have different number of sequences;
#   if so, then check if $ASV_fasta is a subset of $fastasize.size.fasta (by IDs, with stripping ;size=.*).
# If the fasta contains IDs that are not in the ASV table, abort.
size_fa_count=$(grep -c "^>" "$output_dir/$fastasize.size.fasta")
asv_fa_count=$(grep -c "^>" "$ASV_fasta")
if [[ "$size_fa_count" -ne "$asv_fa_count" ]]; then
    printf "\nNote: '%s' has %s sequence(s); ASV table has %s. Verifying that the fasta IDs are a subset of the table.\n" \
        "$ASV_fasta" "$asv_fa_count" "$size_fa_count"
    awk '/^>/ { id = substr($0, 2); sub(/[[:space:];].*/, "", id); print id }' \
        "$output_dir/$fastasize.size.fasta" | sort -u > tempdir/size_fasta_ids.txt
    awk '/^>/ { id = substr($0, 2); sub(/[[:space:];].*/, "", id); print id }' \
        "$ASV_fasta" | sort -u > tempdir/asv_fasta_ids.txt
    missing_ids=$(comm -23 tempdir/asv_fasta_ids.txt tempdir/size_fasta_ids.txt)
    if [[ -n "$missing_ids" ]]; then
        missing_count=$(printf '%s\n' "$missing_ids" | wc -l)
        examples=$(printf '%s\n' "$missing_ids" | head -n 10 | sed 's/^/      - /')
        printf "\nERROR]: '%s' contains %s ID(s) that are not present in the ASV table ('%s').
    The ASV fasta must be a subset of the ASV table (matched by sequence ID, ignoring ';size=...').
    Examples of missing IDs (up to 10):
%s
    >>>Quitting\n" \
            "$(basename "$ASV_fasta")" "$missing_count" "$(basename "$ASV_table_orig")" "$examples" >&2
        end_process
    fi
    printf "OK: all %s fasta ID(s) are present in the ASV table.\n" "$asv_fa_count"
    # Subset the size-annotated fasta to keep only the ASVs that are also in $ASV_fasta,
    # so downstream clustering operates on the same set of sequences as the input fasta.
    printf "Subsetting '%s' to retain only the %s ID(s) present in '%s' ...\n" \
        "$output_dir/$fastasize.size.fasta" "$asv_fa_count" "$ASV_fasta"
    awk 'NR==FNR { ids[$0]=1; next }
         /^>/ { h = substr($0, 2); sub(/[[:space:];].*/, "", h); keep = (h in ids) }
         keep' \
        tempdir/asv_fasta_ids.txt "$output_dir/$fastasize.size.fasta" \
        > tempdir/size_fasta_subset.fasta
    mv tempdir/size_fasta_subset.fasta "$output_dir/$fastasize.size.fasta"
    subset_count=$(grep -c "^>" "$output_dir/$fastasize.size.fasta")
    if [[ "$subset_count" -ne "$asv_fa_count" ]]; then
        printf "\n[ERROR]: After subsetting, '%s' contains %s sequence(s) but expected %s.\n    >>>Quitting\n" \
            "$output_dir/$fastasize.size.fasta" "$subset_count" "$asv_fa_count" >&2
        end_process
    fi
    printf "Done. Size-annotated fasta now contains %s sequence(s).\n" "$subset_count"
fi

# export for R
ASV_fasta_size=$"$output_dir/$fastasize.size.fasta"
export ASV_fasta_size

### Clustering
checkerror=$(vsearch $seqsort \
                    $ASV_fasta_size \
                    $id \
                    $simtype \
                    $strands \
                    $mask \
                    $centroid_in \
                    $maxaccepts \
                    $cores \
                    $otutype $output_dir/OTUs.temp.fasta \
                    --uc $output_dir/OTUs.uc \
                    --fasta_width 0 \
                    --sizein --sizeout 2>&1)
check_app_error


### OTU table creation
printf "Making OTU table ... \n "
Rlog=$(Rscript /scripts/submodules/ASVs2OTUs.R 2>&1)
echo $Rlog > $output_dir/ASVs2OTUs.log 
wait

### Remove size annotations from the OTUs.uc file
sed -i -E 's/;size=[0-9]+//g' "$output_dir/OTUs.uc"

### Discard singleton OTUs
if [[ $remove_singletons == "TRUE"  ]]; then
    printf "Discarding singletons ... \n"
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 2 \
    --sizein --xsize --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error

    mv $output_dir/OTUs.temp.fasta tempdir/
else
    checkerror=$(vsearch \
    --sortbysize $output_dir/OTUs.temp.fasta \
    --minsize 1 \
    --sizein --xsize --fasta_width 0 \
    --output $output_dir/OTUs.fasta 2>&1)
    check_app_error

    mv $output_dir/OTUs.temp.fasta tempdir/
fi

#################################################
### COMPILE FINAL STATISTICS AND README FILES ###
#################################################
printf "\nCleaning up and compiling final stats files ... \n"

#Delete tempdirs if debugger != true
if [[ $debugger != "true" ]]; then
    if [[ -d tempdir ]]; then
        rm -rf tempdir
    fi
    if [[ -d tempdir2 ]]; then
        rm -rf tempdir2
    fi
    if [[ -f $output_dir/ASVs2OTUs.log ]]; then
        rm $output_dir/ASVs2OTUs.log 
    fi
fi

# count features and sequences; outputs variables feature_count, nSeqs, nSample
count_features "$output_dir/OTU_table.txt"
input_ASV_count=$(grep -c "^>" $ASV_fasta_size)
ASV_table_basename=$(basename "$ASV_table")

# Remove temporary cleaned table if it was created
if [[ -n "${ASV_table_nosize}" ]] && [[ -f "${ASV_table_nosize}" ]] && [[ "${debugger:-false}" != "true" ]]; then
    rm -f "$ASV_table_nosize"
fi

end=$(date +%s) 
runtime=$((end-start))

#Make README.txt file
printf "# ASVs clustered to OTUs with vsearch (see 'Core command' below for the used settings).

Start time: $start_time
End time: $(date)
Runtime: $runtime seconds

Clustering formed $feature_count OTUs 
  [input contained $input_ASV_count ASVs (sequences)].

Files in 'ASVs2OTUs' directory:
-------------------------------
# OTUs.fasta            = FASTA formated representative OTU sequences. 
# OTU_table.txt         = OTU distribution table per sample (tab delimited file).
# OTUs.uc               = uclust-like formatted clustering results for OTUs.
# $fastasize.size.fasta = size annotated input sequences [size annotation based on the input table $ASV_table_basename] 

Number of OTUs                       = $feature_count
Number of sequences in the OTU table = $nSeqs
Number of samples in the OTU table   = $nSample

Core command -> 
vsearch $seqsort input.fasta $id $simtype $strands $mask $centroid_in $maxaccepts $cores $otutype OTUs.fasta --fasta_width 0 --sizein
[before clustering, ASV size annotations were fetched from the ASV table]

Total run time was $runtime sec.\n\n
##############################################
###Third-party applications for this process:
#vsearch (version $vsearch_version)
    #citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016) VSEARCH: a versatile open source tool for metagenomics PeerJ 4:e2584
    #https://github.com/torognes/vsearch
################################################" > $output_dir/README.txt

#Done
printf "\nDONE "
printf "Total time: $runtime sec.\n "

#variables for all services
echo "#variables for all services: "
echo "workingDir=$output_dir"
echo "fileFormat=$extension"
echo "readType=single_end"

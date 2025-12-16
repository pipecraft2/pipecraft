#!/bin/bash

# Script to organize pipeline outputs into structured folders
# Usage: bash organize_results.sh <run_id>

RUN_ID="${1:-nano}"
RESULTS_DIR="${RUN_ID}_results"
WORK_DIR="work"

echo "=================================="
echo "Organizing results for: $RUN_ID"
echo "=================================="

# Create organized directory structure
mkdir -p "${RESULTS_DIR}/01_quality_reports"
mkdir -p "${RESULTS_DIR}/02_filtered_sequences"
mkdir -p "${RESULTS_DIR}/03_clusters"
mkdir -p "${RESULTS_DIR}/04_polished_sequences"
mkdir -p "${RESULTS_DIR}/05_its_extracted"
mkdir -p "${RESULTS_DIR}/06_blast_results"
mkdir -p "${RESULTS_DIR}/07_json_results"
mkdir -p "${RESULTS_DIR}/08_processing_logs"

echo ""
echo "Step 1: Moving NanoPlot reports..."
# NanoPlot reports are already in the results directory
if [ -d "${RESULTS_DIR}" ]; then
    for nanoplot_dir in ${RESULTS_DIR}/*_NanoPlot_results; do
        if [ -d "$nanoplot_dir" ]; then
            basename=$(basename "$nanoplot_dir")
            if [ ! -d "${RESULTS_DIR}/01_quality_reports/$basename" ]; then
                mv "$nanoplot_dir" "${RESULTS_DIR}/01_quality_reports/" 2>/dev/null
            fi
        fi
    done
fi
echo "✓ NanoPlot reports organized"

echo ""
echo "Step 2: Collecting filtered sequences..."
find $WORK_DIR -name "*.chopper.fasta.gz" -type f | while read file; do
    sample=$(basename "$file" | sed 's/combined\.\(.*\)\.chopper\.fasta\.gz/\1/')
    cp "$file" "${RESULTS_DIR}/02_filtered_sequences/${sample}.chopper.fasta.gz" 2>/dev/null
done
echo "✓ Filtered sequences collected: $(ls ${RESULTS_DIR}/02_filtered_sequences/*.gz 2>/dev/null | wc -l) files"

echo ""
echo "Step 3: Collecting cluster centroids..."
find $WORK_DIR -name "*.chopper.centeroids.fasta.gz" -type f | while read file; do
    sample=$(basename "$file" | sed 's/combined\.\(.*\)\.chopper\.centeroids\.fasta\.gz/\1/')
    cp "$file" "${RESULTS_DIR}/03_clusters/${sample}.centroids.fasta.gz" 2>/dev/null
done
echo "✓ Cluster centroids collected: $(ls ${RESULTS_DIR}/03_clusters/*.gz 2>/dev/null | wc -l) files"

echo ""
echo "Step 4: Collecting polished sequences (Racon + Medaka)..."
find $WORK_DIR -name "*.racon.fasta" -type f | while read file; do
    sample=$(basename "$file" | sed 's/combined\.\(.*\)\.racon\.fasta/\1/')
    cp "$file" "${RESULTS_DIR}/04_polished_sequences/${sample}.racon.fasta" 2>/dev/null
done

find $WORK_DIR -name "consensus.fasta" -type f | while read file; do
    medaka_dir=$(dirname "$file")
    sample=$(basename "$medaka_dir" | sed 's/\(.*\)_medaka_output/\1/')
    if [ ! -z "$sample" ]; then
        cp "$file" "${RESULTS_DIR}/04_polished_sequences/${sample}.medaka.consensus.fasta" 2>/dev/null
    fi
done
echo "✓ Polished sequences collected: $(ls ${RESULTS_DIR}/04_polished_sequences/*.fasta 2>/dev/null | wc -l) files"

echo ""
echo "Step 5: Collecting ITS extracted sequences..."
find $WORK_DIR -name "*.after_itsx.fasta" -type f | while read file; do
    sample=$(basename "$file" .after_itsx.fasta)
    cp "$file" "${RESULTS_DIR}/05_its_extracted/${sample}.its.fasta" 2>/dev/null
done
echo "✓ ITS sequences collected: $(ls ${RESULTS_DIR}/05_its_extracted/*.fasta 2>/dev/null | wc -l) files"

echo ""
echo "Step 6: Collecting BLAST results..."
find $WORK_DIR -name "*_blastn_results.tsv" -type f | while read file; do
    sample=$(basename "$file" _blastn_results.tsv)
    cp "$file" "${RESULTS_DIR}/06_blast_results/${sample}.blast.tsv" 2>/dev/null
done
echo "✓ BLAST results collected: $(ls ${RESULTS_DIR}/06_blast_results/*.tsv 2>/dev/null | wc -l) files"

echo ""
echo "Step 7: Collecting JSON results..."
find $WORK_DIR -name "*.results.json" -type f | while read file; do
    sample=$(basename "$file" .results.json)
    cp "$file" "${RESULTS_DIR}/07_json_results/${sample}.json" 2>/dev/null
done
echo "✓ JSON results collected: $(ls ${RESULTS_DIR}/07_json_results/*.json 2>/dev/null | wc -l) files"

echo ""
echo "Step 8: Collecting processing logs..."
find $WORK_DIR -name "processing.log" -type f | while read file; do
    processing_dir=$(dirname "$file")
    sample=$(basename "$processing_dir" | sed 's/processing_//')
    if [ ! -z "$sample" ]; then
        cp "$file" "${RESULTS_DIR}/08_processing_logs/${sample}.processing.log" 2>/dev/null
    fi
done
echo "✓ Processing logs collected: $(ls ${RESULTS_DIR}/08_processing_logs/*.log 2>/dev/null | wc -l) files"

# Create a summary README
cat > "${RESULTS_DIR}/README.txt" << EOF
====================================================================
FunBarONT Analysis Results - Run ID: $RUN_ID
====================================================================
Generated on: $(date)

DIRECTORY STRUCTURE:
--------------------

${RUN_ID}.results.xlsx
    → Final aggregated results with taxonomic assignments

01_quality_reports/
    → NanoPlot quality assessment reports for each sample
    → Open the HTML files in a web browser

02_filtered_sequences/
    → Quality-filtered sequences (Chopper)
    → Format: FASTA.GZ

03_clusters/
    → Clustered sequence centroids (VSEARCH)
    → Format: FASTA.GZ

04_polished_sequences/
    → Racon polished sequences (.racon.fasta)
    → Medaka polished consensus sequences (.medaka.consensus.fasta)
    → Format: FASTA

05_its_extracted/
    → ITS regions extracted using ITSx
    → Format: FASTA

06_blast_results/
    → BLAST results against UNITE database
    → Format: TSV (tab-separated values)
    → Columns: qseqid, sseqid, pident, qcovs, evalue, qlen, slen

07_json_results/
    → Detailed JSON results for each sample
    → Contains cluster data, BLAST hits, and metadata

08_processing_logs/
    → Processing logs for each sample
    → Contains timestamps and step-by-step execution details

====================================================================
QUICK START:
--------------------
1. Open ${RUN_ID}.results.xlsx for final taxonomic assignments
2. View quality reports in 01_quality_reports/
3. Check individual sample details in 07_json_results/

====================================================================
EOF

echo ""
echo "=================================="
echo "✓ Organization complete!"
echo "=================================="
echo ""
echo "Results organized in: $RESULTS_DIR/"
echo ""
echo "Directory summary:"
echo "  - Final Excel results: ${RUN_ID}.results.xlsx"
echo "  - Quality reports: 01_quality_reports/ ($(ls ${RESULTS_DIR}/01_quality_reports/ 2>/dev/null | wc -l) samples)"
echo "  - Filtered sequences: 02_filtered_sequences/ ($(ls ${RESULTS_DIR}/02_filtered_sequences/ 2>/dev/null | wc -l) files)"
echo "  - Clusters: 03_clusters/ ($(ls ${RESULTS_DIR}/03_clusters/ 2>/dev/null | wc -l) files)"
echo "  - Polished sequences: 04_polished_sequences/ ($(ls ${RESULTS_DIR}/04_polished_sequences/ 2>/dev/null | wc -l) files)"
echo "  - ITS extracted: 05_its_extracted/ ($(ls ${RESULTS_DIR}/05_its_extracted/ 2>/dev/null | wc -l) files)"
echo "  - BLAST results: 06_blast_results/ ($(ls ${RESULTS_DIR}/06_blast_results/ 2>/dev/null | wc -l) files)"
echo "  - JSON results: 07_json_results/ ($(ls ${RESULTS_DIR}/07_json_results/ 2>/dev/null | wc -l) files)"
echo "  - Processing logs: 08_processing_logs/ ($(ls ${RESULTS_DIR}/08_processing_logs/ 2>/dev/null | wc -l) files)"
echo ""
echo "Read README.txt in the results directory for more information."
echo "=================================="

#!/bin/bash

############################
# FunBarONT Pipeline Wrapper
############################
# Oxford Nanopore Technologies fungal barcoding pipeline
# Processes ONT basecaller output into high-quality ITS sequences

# Check if required tools are available
command -v nextflow >/dev/null 2>&1 || { echo "ERROR: nextflow not found"; exit 1; }
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found"; exit 1; }

echo "=========================================="
echo "FunBarONT Pipeline Starting"
echo "=========================================="
echo "Date: $(date)"
echo ""

# Validate input directory
if [ ! -d "/Input" ]; then
    echo "ERROR: No input directory mounted at /Input"
    exit 1
fi

if [ ! -d "/sequences" ]; then
    echo "ERROR: No sequences directory mounted at /sequences"
    exit 1
fi

# Read configuration from JSON file
CONFIG_FILE="/scripts/FunBarONTConfig.json"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

echo "Reading configuration from: $CONFIG_FILE"

DATABASE_FILE=$(jq -r '.database_file' "$CONFIG_FILE")
BLASTDB_PATH=$(jq -r '.blastdb_path' "$CONFIG_FILE")
RUN_ID=$(jq -r '.run_id' "$CONFIG_FILE")
MEDAKA_MODEL=$(jq -r '.medaka_model // "r1041_e82_400bps_hac_variant_v4.3.0"' "$CONFIG_FILE")
USE_ITSX=$(jq -r '.use_itsx // 1' "$CONFIG_FILE")
OUTPUT_ALL_POLISHED=$(jq -r '.output_all_polished_seqs // 0' "$CONFIG_FILE")
REL_ABU_THRESHOLD=$(jq -r '.rel_abu_threshold // 10' "$CONFIG_FILE")
CPU_THREADS=$(jq -r '.cpu_threads // 8' "$CONFIG_FILE")
CHOPPER_QUALITY=$(jq -r '.chopper_quality // 10' "$CONFIG_FILE")
CHOPPER_MIN_LENGTH=$(jq -r '.chopper_min_read_length // 150' "$CONFIG_FILE")
CHOPPER_MAX_LENGTH=$(jq -r '.chopper_max_read_length // 1000' "$CONFIG_FILE")
VSEARCH_CLUSTER_ID=$(jq -r '.vsearch_cluster_id // 0.95' "$CONFIG_FILE")

echo ""
echo "Configuration:"
echo "  Database File: $DATABASE_FILE"
echo "  BLAST DB Path: $BLASTDB_PATH"
echo "  Run ID: $RUN_ID"
echo "  Medaka Model: $MEDAKA_MODEL"
echo "  Use ITSx: $USE_ITSX"
echo "  Output All Polished: $OUTPUT_ALL_POLISHED"
echo "  Rel Abu Threshold: $REL_ABU_THRESHOLD"
echo "  CPU Threads: $CPU_THREADS"
echo "  Chopper Min/Max Length: $CHOPPER_MIN_LENGTH / $CHOPPER_MAX_LENGTH"
echo "  VSEARCH Cluster ID: $VSEARCH_CLUSTER_ID"
echo ""

# Check if input directory exists and has data
if [ ! -d "/Input" ]; then
    echo "ERROR: Input directory not found: /Input"
    exit 1
fi

# Count input files
FASTQ_COUNT=$(find /Input -maxdepth 1 -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)
echo "Found $FASTQ_COUNT FASTQ files in input directory"

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in /Input"
    exit 1
fi

# Create BLAST database directory
mkdir -p "$BLASTDB_PATH"
mkdir -p "$BLASTDB_PATH/unite"

# Check if database file exists
if [ ! -f "$DATABASE_FILE" ]; then
    echo "ERROR: Database file not found: $DATABASE_FILE"
    exit 1
fi

echo ""
echo "=========================================="
echo "Preparing BLAST Database"
echo "=========================================="
echo ""

# Check if BLAST database already exists
if [ -f "$BLASTDB_PATH/unite/unite.nhr" ]; then
    echo "BLAST database already exists. Skipping makeblastdb..."
else
    echo "Creating BLAST database from: $DATABASE_FILE"
    
    # Copy database file to BLAST directory
    cp "$DATABASE_FILE" "$BLASTDB_PATH/unite/unite.fasta"
    
    # Create BLAST database (without -parse_seqids to allow long headers in UNITE)
    makeblastdb \
        -in "$BLASTDB_PATH/unite/unite.fasta" \
        -dbtype nucl \
        -out "$BLASTDB_PATH/unite/unite" \
        -blastdb_version 5
    
    MAKEBLASTDB_EXIT=$?
    
    if [ $MAKEBLASTDB_EXIT -ne 0 ]; then
        echo "ERROR: makeblastdb failed with exit code $MAKEBLASTDB_EXIT"
        exit $MAKEBLASTDB_EXIT
    fi
    
    echo "BLAST database created successfully"
fi

echo ""
WORK_DIR="/sequences/work"
RESULTS_DIR="/sequences/${RUN_ID}_results"
mkdir -p "$WORK_DIR"
mkdir -p "$RESULTS_DIR"

# Set Nextflow home to writable location
export NXF_HOME="/sequences/.nextflow"
mkdir -p "$NXF_HOME"

echo ""
echo "Working directory: $WORK_DIR"
echo "Results directory: $RESULTS_DIR"
echo ""

# Start time tracking
start_time=$(date +%s)

echo "=========================================="
echo "Running Nextflow Pipeline"
echo "=========================================="
echo ""

# Run Nextflow pipeline from /sequences directory (writable) instead of /scripts
# This avoids permission issues with .nextflow directory creation
cd /sequences || exit 1

# Run Nextflow pipeline
nextflow run /scripts/FunBarONT/main.nf \
  --FASTQ_DIRECTORY "/Input" \
  --BLASTDB_PATH "$BLASTDB_PATH" \
  --RUN_ID "$RUN_ID" \
  --MEDAKA_MODEL "$MEDAKA_MODEL" \
  --USE_ITSX "$USE_ITSX" \
  --OUTPUT_ALL_POLISHED_SEQS "$OUTPUT_ALL_POLISHED" \
  --CHOPPER_MIN_READ_LENGTH "$CHOPPER_MIN_LENGTH" \
  --CHOPPER_MAX_READ_LENGTH "$CHOPPER_MAX_LENGTH" \
  --REL_ABU_THRESHOLD "$REL_ABU_THRESHOLD" \
  --CPU_THREADS "$CPU_THREADS" \
  -work-dir "$WORK_DIR" \
  -resume

NEXTFLOW_EXIT=$?

echo ""
echo "=========================================="
echo "Nextflow Pipeline Completed"
echo "=========================================="
echo "Exit code: $NEXTFLOW_EXIT"
echo ""

if [ $NEXTFLOW_EXIT -ne 0 ]; then
    echo "ERROR: Nextflow pipeline failed with exit code $NEXTFLOW_EXIT"
    echo "Check the Nextflow logs for details"
    exit $NEXTFLOW_EXIT
fi

# Calculate runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
runtime_min=$((runtime / 60))
runtime_sec=$((runtime % 60))

echo "Pipeline runtime: ${runtime_min}m ${runtime_sec}s"
echo ""

# Generate summary README
echo "=========================================="
echo "Generating Pipeline Summary"
echo "=========================================="

README_FILE="$RESULTS_DIR/README.md"

cat > "$README_FILE" << EOF
# FunBarONT Pipeline Results

## Run Information

- **Run ID**: $RUN_ID
- **Date**: $(date)
- **Runtime**: ${runtime_min}m ${runtime_sec}s
- **Pipeline**: FunBarONT - Oxford Nanopore Technologies Fungal Barcoding

## Input Data

- **Input Directory**: /Input
- **Input FASTQ files**: $FASTQ_COUNT
- **Database File**: $DATABASE_FILE
- **BLAST Database**: $BLASTDB_PATH

## Pipeline Parameters

- **Medaka Model**: $MEDAKA_MODEL
- **ITSx Extraction**: $([ "$USE_ITSX" = "1" ] && echo "Enabled" || echo "Disabled")
- **Output All Polished**: $([ "$OUTPUT_ALL_POLISHED" = "1" ] && echo "Enabled" || echo "Disabled")
- **Rel Abu Threshold**: $REL_ABU_THRESHOLD%
- **CPU Threads**: $CPU_THREADS
- **Chopper Quality**: $CHOPPER_QUALITY
- **Chopper Min Length**: $CHOPPER_MIN_LENGTH bp
- **Chopper Max Length**: $CHOPPER_MAX_LENGTH bp
- **VSEARCH Cluster ID**: $VSEARCH_CLUSTER_ID

## Pipeline Steps

### 1. Quality Control (NanoPlot)
- **Tool**: NanoPlot
- **Purpose**: Assess quality metrics of ONT sequencing data
- **Output**: Quality plots and statistics for each barcode

### 2. Quality Filtering (chopper)
- **Tool**: chopper
- **Purpose**: Filter reads based on quality scores
- **Parameters**: 
  - Min quality = $CHOPPER_QUALITY
  - Min length = $CHOPPER_MIN_LENGTH bp
  - Max length = $CHOPPER_MAX_LENGTH bp
- **Output**: Filtered FASTQ files

### 3. Read Mapping (minimap2)
- **Tool**: minimap2
- **Purpose**: Map reads to reference sequences for clustering
- **Output**: BAM alignment files

### 4. Consensus Calling (racon + medaka)
- **Tool**: racon (1st pass), medaka (polishing)
- **Purpose**: Generate high-accuracy consensus sequences
- **Medaka Model**: $MEDAKA_MODEL
- **Output**: Polished consensus FASTA

### 5. Clustering (vsearch)
- **Tool**: vsearch
- **Purpose**: Cluster similar sequences
- **Parameters**: Cluster identity = $VSEARCH_CLUSTER_ID
- **Output**: Representative sequences (OTUs)

EOF

if [ "$USE_ITSX" = "1" ]; then
cat >> "$README_FILE" << EOF
### 6. ITS Extraction (ITSx)
- **Tool**: ITSx
- **Purpose**: Extract ITS1, 5.8S, and ITS2 regions from fungal sequences
- **Output**: ITS region FASTA files

EOF
fi

cat >> "$README_FILE" << EOF
### $( [ "$USE_ITSX" = "1" ] && echo "7" || echo "6" ). Taxonomy Assignment (BLAST)
- **Tool**: BLAST+ (blastn)
- **Database**: $BLASTDB_PATH
- **Purpose**: Assign taxonomy to consensus sequences
- **Output**: BLAST results with top hits

## Output Files

### Main Results
- \`results_table.xlsx\` - Excel table with all results (sequence info, taxonomy, quality metrics)
- \`consensus_sequences.fasta\` - All consensus sequences
- \`taxonomy_assignments.txt\` - Detailed taxonomy assignments

### Quality Control
- \`nanoplot/\` - Quality assessment plots and statistics per barcode

### Intermediate Files
- \`filtered_reads/\` - Quality-filtered FASTQ files
- \`alignments/\` - Minimap2 alignment files
- \`consensus/\` - Consensus sequences per barcode
- \`clusters/\` - Clustering results

EOF

if [ "$USE_ITSX" = "1" ]; then
cat >> "$README_FILE" << EOF
### ITS Extraction
- \`itsx/\` - ITS region sequences (ITS1, 5.8S, ITS2)

EOF
fi

cat >> "$README_FILE" << EOF
### BLAST Results
- \`blast/\` - BLAST output files
- \`taxonomy/\` - Parsed taxonomy assignments

## Sequence Statistics

EOF

# Try to count sequences in output files
if [ -f "$RESULTS_DIR/consensus_sequences.fasta" ]; then
    CONSENSUS_COUNT=$(grep -c "^>" "$RESULTS_DIR/consensus_sequences.fasta" || echo "0")
    echo "- **Total consensus sequences**: $CONSENSUS_COUNT" >> "$README_FILE"
fi

cat >> "$README_FILE" << EOF

## Tools and Citations

### Core Pipeline Tools

**FunBarONT**
- Citation: [Provide FunBarONT citation]
- GitHub: https://github.com/mdziurzynski/ont_fungal_barcoding_pipeline

**Nextflow**
- Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol 35:316-319
- Website: https://www.nextflow.io/

**NanoPlot**
- De Coster W, et al. (2018) NanoPack: visualizing and processing long-read sequencing data. Bioinformatics 34:2666-2669
- GitHub: https://github.com/wdecoster/NanoPlot

**chopper**
- De Coster W, Rademakers R (2023) NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics 39:btad311
- GitHub: https://github.com/wdecoster/chopper

**minimap2**
- Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100
- GitHub: https://github.com/lh3/minimap2

**racon**
- Vaser R, et al. (2017) Fast and accurate de novo genome assembly from long uncorrected reads. Genome Res 27:737-746
- GitHub: https://github.com/isovic/racon

**medaka**
- Oxford Nanopore Technologies
- GitHub: https://github.com/nanoporetech/medaka

**vsearch**
- Rognes T, et al. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584
- GitHub: https://github.com/torognes/vsearch

**ITSx**
- Bengtsson-Palme J, et al. (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods Ecol Evol 4:914-919

**BLAST+**
- Camacho C, et al. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421
- Website: https://blast.ncbi.nlm.nih.gov/

### Analysis Platform

**PipeCraft**
- Anslan S, et al. (2017) PipeCraft: flexible open-source toolkit for bioinformatics analysis of custom high-throughput amplicon sequencing data. Mol Ecol Resour 17:e234-e240
- GitHub: https://github.com/pipecraft2/pipecraft

---

**Pipeline completed successfully**

For questions or issues, please refer to:
- FunBarONT documentation: https://github.com/mdziurzynski/ont_fungal_barcoding_pipeline
- PipeCraft documentation: https://pipecraft2-manual.readthedocs.io/

EOF

echo "Summary README generated: $README_FILE"
echo ""

# Copy results to sequences directory
echo "=========================================="
echo "Copying Results"
echo "=========================================="

if [ -d "results" ]; then
    echo "Copying results from Nextflow output..."
    cp -r results/* "$RESULTS_DIR/" 2>/dev/null || echo "Note: Some result files may not exist yet"
fi

# Set ownership if HOST_UID and HOST_GID are provided
if [ ! -z "$HOST_UID" ] && [ ! -z "$HOST_GID" ]; then
    echo "Setting ownership of results to $HOST_UID:$HOST_GID"
    chown -R "$HOST_UID:$HOST_GID" /sequences
fi

echo ""
echo "=========================================="
echo "FunBarONT Pipeline Complete!"
echo "=========================================="
echo ""
echo "Results location: $RESULTS_DIR"
echo "README: $README_FILE"
echo ""
echo "Pipeline finished successfully at $(date)"
echo "=========================================="

exit 0

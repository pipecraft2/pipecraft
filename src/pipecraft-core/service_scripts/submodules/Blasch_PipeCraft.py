#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import shutil
import subprocess
import time
import csv
import multiprocessing
import psutil
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import defaultdict

###############################################################################
#               False positive chimera detection and recovery                 #     
#               for metabarcoding and environmental DNA (eDNA)                #
#                                                                             #
# Description:                                                                #
#   This script processes BLAST XML results to identify, classify, and        #
#   recover false positive chimeric sequences from metabarcoding or eDNA      #
#   datasets. Sequences are categorized into three groups: non-chimeric,      #
#   absolute chimeras, and borderline sequences, based on identity and        #
#   coverage thresholds.                                                      #
#                                                                             #
# Author:        ALI HAKIMZADEH                                               #
# Version:       1.0                                                          #
# Date:          2025-04-01                                                   #
###############################################################################

###############################################################################
#                                Logging Setup                                #
###############################################################################

logger = logging.getLogger("chimera_recovery")
logger.setLevel(logging.INFO)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

###############################################################################
#                               Helper Functions                              #
###############################################################################

def log_system_usage():
    """Logs CPU and RAM usage."""
    process = psutil.Process(os.getpid())
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_info = process.memory_info()
    logger.info(f"CPU usage: {cpu_percent}%")
    logger.info(f"Memory usage: {memory_info.rss / (1024 ** 2):.2f} MB")


def clean_directory(dir_path):
    """Remove all contents in a directory without deleting the directory itself."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)
        return
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logger.error(f"Failed to delete {file_path}. Reason: {e}")


def check_blast_db_exists(db_prefix):
    """Check if a BLAST DB with prefix db_prefix has .nin/.nsq/.nhr files."""
    for ext in [".nin", ".nsq", ".nhr"]:
        if os.path.isfile(db_prefix + ext):
            return True
    return False


def create_blast_db_if_needed(fasta_or_db_prefix, db_dir, db_name):
    """
    Check if 'fasta_or_db_prefix' is already a valid BLAST DB.
    If not, assume it's a FASTA and try to create a DB.
    Return the final DB prefix path if successful, or exit on failure.
    """
    if check_blast_db_exists(fasta_or_db_prefix):
        # It's already a valid DB
        return fasta_or_db_prefix

    # Otherwise, treat it as FASTA and build a new DB
    logger.info(f"'{fasta_or_db_prefix}' doesn't look like a valid BLAST DB. Trying to create DB...")

    os.makedirs(db_dir, exist_ok=True)
    out_prefix = os.path.join(db_dir, db_name)

    # Validate the FASTA file with BioPython
    if not validate_fasta_with_biopython(fasta_or_db_prefix):
        logger.error(f"Reference file '{fasta_or_db_prefix}' is not a valid FASTA.")
        sys.exit(1)

    cmd = [
        "makeblastdb",
        "-in", fasta_or_db_prefix,
        "-dbtype", "nucl",
        "-out", out_prefix
    ]
    logger.info(f"Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create BLAST DB from '{fasta_or_db_prefix}': {e}")
        sys.exit(1)

    if check_blast_db_exists(out_prefix):
        logger.info(f"Created BLAST DB at {out_prefix}")
        return out_prefix
    else:
        logger.error(f"Failed to create a valid DB at {out_prefix}. Exiting.")
        sys.exit(1)


def validate_fasta_with_biopython(fasta_path):
    """
    Validate a FASTA file using BioPython by attempting to parse it.
    Returns True if parsing succeeds and at least one record is found,
    otherwise False.
    """
    if not os.path.isfile(fasta_path):
        logger.error(f"File not found: {fasta_path}")
        return False

    try:
        count = 0
        with open(fasta_path, "r") as handle:
            for _ in SeqIO.parse(handle, "fasta"):
                count += 1
        if count == 0:
            logger.error(f"No sequences found in FASTA: {fasta_path}")
            return False
        return True
    except Exception as e:
        logger.error(f"Error parsing FASTA {fasta_path}: {e}")
        return False

###############################################################################
#                      STEP 1: Create Self-Databases                          #
###############################################################################

def create_self_databases(self_fasta_dir, self_db_dir):
    """
    Creates self-databases from FASTA files in self_fasta_dir.
    Auto-detects FASTA files and creates databases for chimera analysis.
    Prioritizes original sample files over .chimeras files.
    Exits if creation fails for any file.
    """
    logger.info("=== Step 1: Creating Self-Databases ===")
    os.makedirs(self_db_dir, exist_ok=True)

    # Look for FASTA files with various extensions
    fasta_extensions = [".fasta", ".fa", ".fas", ".fna"]
    all_files = os.listdir(self_fasta_dir)
    
    # Separate sample files from chimeras files
    sample_files = []
    chimeras_files = []
    
    for file in all_files:
        if any(file.endswith(ext) for ext in fasta_extensions):
            if '.chimeras.' in file:
                chimeras_files.append(file)
            else:
                sample_files.append(file)
    
    # Prioritize sample files over chimeras files
    if sample_files:
        fasta_files = sample_files
        logger.info(f"Found {len(fasta_files)} sample FASTA files for database creation.")
    elif chimeras_files:
        fasta_files = chimeras_files
        logger.info(f"No sample files found. Using {len(fasta_files)} .chimeras files for database creation.")
    else:
        logger.error(f"No FASTA files found in {self_fasta_dir} with extensions {fasta_extensions}.")
        logger.error("Please ensure that sample FASTA files or .chimeras files are present in the working directory.")
        sys.exit(1)

    logger.info(f"Creating self-databases from {len(fasta_files)} FASTA files...")

    for fasta_file in fasta_files:
        full_path = os.path.join(self_fasta_dir, fasta_file)

        # Validate with BioPython
        if not validate_fasta_with_biopython(full_path):
            logger.error(f"Invalid FASTA file: {full_path}")
            sys.exit(1)

        # Extract base name (remove extension and any .chimeras suffix)
        base_name = os.path.splitext(fasta_file)[0]
        if base_name.endswith('.chimeras'):
            base_name = base_name[:-9]  # Remove .chimeras suffix
            
        out_dir = os.path.join(self_db_dir, base_name)
        os.makedirs(out_dir, exist_ok=True)

        db_prefix = os.path.join(out_dir, base_name)
        
        # Check if database already exists
        if check_blast_db_exists(db_prefix):
            logger.info(f"Self BLAST database already exists for {base_name}, skipping creation.")
            continue
            
        cmd = [
            "makeblastdb",
            "-in", full_path,
            "-dbtype", "nucl",
            "-out", db_prefix
        ]

        logger.info(f"Creating self BLAST database for {fasta_file} -> {base_name} ...")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating database for {fasta_file}: {e}")
            logger.error(f"makeblastdb stderr: {e.stderr if hasattr(e, 'stderr') else 'No stderr'}")
            sys.exit(1)  # exit on failure

        # Verify database was created successfully
        if check_blast_db_exists(db_prefix):
            logger.info(f"Successfully created self DB for {base_name}")
        else:
            logger.error(f"Failed to create valid database for {base_name}")
            sys.exit(1)

    logger.info("All self-databases have been created.\n")


###############################################################################
#                          STEP 2: Reference DB Setup                         #
###############################################################################

def handle_reference_db(reference_db_arg, output_dir):
    """
    If reference_db_arg is empty, return empty string (no reference DB).
    If reference_db_arg is not empty:
      - If it is already a valid DB, return it.
      - Else try to create a DB from it (assuming it's FASTA),
        storing it in output_dir/reference_db/reference
    """
    if not reference_db_arg:
        logger.info("No reference DB provided.")
        return ""
    
    logger.info(f"Reference database argument provided: {reference_db_arg}")
    
    # Check if reference_db_arg is actually a file that exists
    if not os.path.exists(reference_db_arg):
        logger.error(f"Reference database path does not exist: {reference_db_arg}")
        logger.error("Please provide a valid file path for the reference database.")
        sys.exit(1)
        
    reference_db_dir = os.path.join(output_dir, "reference_db")
    db_name = "reference"

    final_db_prefix = create_blast_db_if_needed(reference_db_arg, reference_db_dir, db_name)
    
    # Extra validation of the final database
    if not check_blast_db_exists(final_db_prefix):
        logger.error(f"Final reference database validation failed: {final_db_prefix}")
        logger.error("Unable to create or validate reference database. Check file permissions and disk space.")
        sys.exit(1)
        
    logger.info(f"Reference database successfully validated: {final_db_prefix}")
    return final_db_prefix


###############################################################################
#                     STEP 3: Run BLAST on .chimeras.fasta                    #
###############################################################################

def run_blast_analysis(
    input_chimeras_dir,
    self_db_dir,
    reference_db_prefix,
    output_dir,
    threads
):
    """
    For each .chimeras.fasta in input_chimeras_dir, run BLAST against:
      - reference_db_prefix (if not empty)
      - corresponding self-database
    Store XML in output_dir/xml/<base_name>_blast_results.xml
    Exit on first BLAST failure.
    """
    logger.info("=== Step 3: Running BLAST Analyses ===")
    
    # Debug reference database information
    if reference_db_prefix:
        logger.info(f"Reference database will be used: {reference_db_prefix}")
        # Check if reference database exists and is readable
        if check_blast_db_exists(reference_db_prefix):
            logger.info(f"Reference DB validation successful: Found BLAST database files for {reference_db_prefix}")
        else:
            logger.warning(f"Reference DB validation FAILED: No BLAST database files found at {reference_db_prefix}")
    else:
        logger.info("No reference database provided, will use only self-databases")
    xml_dir = os.path.join(output_dir, "xml")
    os.makedirs(xml_dir, exist_ok=True)

    # Collect .chimeras files with various extensions
    fasta_extensions = [".fasta", ".fa", ".fas", ".fna"]
    chimeric_files = []
    
    for ext in fasta_extensions:
        chimeric_files.extend([f for f in os.listdir(input_chimeras_dir) if f.endswith(f".chimeras{ext}")])
    
    if not chimeric_files:
        logger.error(f"No .chimeras files found in {input_chimeras_dir} with extensions {fasta_extensions}.")
        logger.error("Please ensure chimera detection has been run prior to BlasCh analysis.")
        sys.exit(1)
    
    logger.info(f"Found {len(chimeric_files)} chimeras files for BLAST analysis.")

    for chim_file in chimeric_files:
        chimera_file = os.path.join(input_chimeras_dir, chim_file)

        # Validate the chimera_file with BioPython
        if not validate_fasta_with_biopython(chimera_file):
            logger.error(f"Invalid chimeras file: {chimera_file}")
            sys.exit(1)

        base_name = os.path.splitext(chim_file)[0]  # e.g., "ERR6454463.chimeras"

        # The original sample name might be base_name.replace(".chimeras", "")
        # e.g. "ERR6454463"
        raw_base = base_name.replace(".chimeras", "")
        self_subdir = os.path.join(self_db_dir, raw_base)
        self_db = os.path.join(self_subdir, raw_base)

        if not check_blast_db_exists(self_db):
            logger.error(f"Self DB not found or invalid for {chim_file}: {self_db}")
            sys.exit(1)

        # Output XML => <base_name>_blast_results.xml
        xml_file = os.path.join(xml_dir, f"{base_name}_blast_results.xml")

        # Build the BLAST command
        cmd = [
            "blastn",
            "-query", chimera_file,
            # Don't include -db yet - we'll add it based on database options
            "-word_size", "7",
            "-task", "blastn",
            "-num_threads", str(threads),
            "-outfmt", "5",
            "-evalue", "0.001",
            "-strand", "both",
            "-max_target_seqs", "10",
            "-max_hsps", "9",
            "-out", xml_file
        ]
        
        # Add database arguments correctly - must be a single -db followed by all databases
        if reference_db_prefix:
            # In BLAST+ when passing multiple databases, they must be specified as one space-separated string
            # This allows BLAST to search in all databases at once
            db_value = f"{reference_db_prefix} {self_db}"
            cmd.extend(["-db", db_value])
            logger.info(f"Running BLAST for {chim_file} against multiple DBs: {db_value}")
            # Debug: Print the exact command being run
            logger.info(f"Full BLAST command will be: {' '.join(cmd)}")
        else:
            cmd.extend(["-db", self_db])
            logger.info(f"Running BLAST for {chim_file} against single DB: {self_db}")
        
        # Debug: log the exact command being run
        logger.info(f"Running BLAST command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST failed for {chim_file}: {e}")
            sys.exit(1)  # exit on BLAST failure

        logger.info(f"Completed BLAST analysis for {base_name}")

    logger.info("All BLAST analyses completed.\n")


###############################################################################
#             STEP 4: Chimera Detection and Sequence Recovery                 #
###############################################################################

def load_fasta_sequences(fasta_file):
    """Load sequences from a FASTA file into a dict."""
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    seqs = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        seqs[rec.id] = str(rec.seq)
    return seqs

def extract_query_id(blast_query_def):
    """Extract the query ID from the BLAST query definition."""
    return blast_query_def

def is_self_hit(hit_def):
    """Check if a hit is from the same sample based on ';size=' in the def line."""
    return ';size=' in hit_def

def extract_taxonomy(hit_def):
    """Extract taxonomy from the first semicolon onward, else 'Unclassified'."""
    try:
        return hit_def.split(';', 1)[1]
    except IndexError:
        return "Unclassified"

# Classification thresholds
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
BORDERLINE_COVERAGE_THRESHOLD = 89.0
MULTIPLE_ALIGNMENT_COVERAGE_LIMIT = 85.0  # Fixed threshold for multiple alignment coverage

def analyze_blast_hits(blast_record, query_id):
    """
    Examine the BLAST hits (alignments).
    If the first non-self alignment has multiple HSPs => forced chimeric.
    """
    hits_info = []
    if not blast_record.alignments:
        return hits_info

    first_alignment = blast_record.alignments[0]
    first_hit_id = extract_query_id(first_alignment.hit_def)

    hit_to_check = first_alignment
    alt_hit_id = first_hit_id

    if first_hit_id == query_id and len(blast_record.alignments) > 1:
        hit_to_check = blast_record.alignments[1]
        alt_hit_id = extract_query_id(hit_to_check.hit_def)

    # Check multiple HSPs in first non-self alignment
    if alt_hit_id != query_id and len(hit_to_check.hsps) > 1:
        # Calculate best coverage among HSPs
        best_hsp = max(hit_to_check.hsps, key=lambda hsp: (hsp.align_length / blast_record.query_length) * 100)
        best_coverage = min((best_hsp.align_length / blast_record.query_length) * 100, 100)
        
        hits_info.append({
            "hit_id": alt_hit_id,
            "identity": 100.0,
            "coverage": best_coverage,
            "is_same_sample": False,
            "taxonomy": extract_taxonomy(hit_to_check.hit_def),
            "force_chimeric": best_coverage <= MULTIPLE_ALIGNMENT_COVERAGE_LIMIT,
            "multiple_hsps": True
        })
        return hits_info

    # Otherwise, gather hits
    db_hit_count = 0
    self_hit_count = 0
    
    for aln in blast_record.alignments:
        if not aln.hsps:
            continue
        candidate_id = extract_query_id(aln.hit_def)
        if candidate_id == query_id:  # self alignment
            continue

        same_sample = is_self_hit(aln.hit_def)
        if same_sample:
            self_hit_count += 1
        else:
            db_hit_count += 1
            
        taxonomy = extract_taxonomy(aln.hit_def) if not same_sample else "Self"

        best_hsp = max(
            aln.hsps,
            key=lambda hsp: (hsp.align_length / blast_record.query_length) * 100
        )
        identity_pct = 0.0
        coverage_pct = 0.0
        if best_hsp.align_length > 0:
            identity_pct = (best_hsp.identities / best_hsp.align_length) * 100
            coverage_pct = min(
                (best_hsp.align_length / blast_record.query_length) * 100, 100
            )

        hits_info.append({
            "hit_id": candidate_id,
            "identity": identity_pct,
            "coverage": coverage_pct,
            "is_same_sample": same_sample,
            "taxonomy": taxonomy,
            "force_chimeric": False,
            "multiple_hsps": False
        })

    # Debug logging
    if db_hit_count > 0 or self_hit_count > 0:
        logger.info(f"Query {query_id}: Found {db_hit_count} reference DB hits, {self_hit_count} self hits")

    return hits_info

def classify_sequence(
    hits_info,
    high_identity_threshold,
    high_coverage_threshold,
    borderline_identity_threshold,
    borderline_coverage_threshold
):
    """
    Apply classification logic:
      1. If forced chimeric => chimeric
      2. If no hits => non_chimeric
      3. If only self hits => chimeric
      4. High quality => non_chimeric
      5. Borderline => non_chimeric
      6. Else borderline or chimeric
    """
    if not hits_info:
        return "non_chimeric", "No significant non-self hits"
    if any(h["force_chimeric"] for h in hits_info):
        return "chimeric", "Multiple HSPs in first non-self hit with coverage <= 85%"

    db_hits = [h for h in hits_info if not h["is_same_sample"]]
    self_hits = [h for h in hits_info if h["is_same_sample"]]

    # Debug logging for classification
    logger.debug(f"Classification input: {len(db_hits)} DB hits, {len(self_hits)} self hits")
    if db_hits:
        logger.debug(f"Best DB hit: identity={db_hits[0]['identity']:.1f}%, coverage={db_hits[0]['coverage']:.1f}%")

    # If no db hits => chimeric
    if not db_hits and self_hits:
        return "chimeric", "Only self-hits, no DB hits"

    if db_hits:
        # High-quality
        high_quality_hits = [
            h for h in db_hits
            if h["identity"] >= high_identity_threshold
               and h["coverage"] >= high_coverage_threshold
        ]
        if high_quality_hits:
            return "non_chimeric", "High-quality DB match"

        # Borderline
        borderline_hits = [
            h for h in db_hits
            if h["identity"] >= borderline_identity_threshold
               and h["coverage"] >= borderline_coverage_threshold
        ]
        if borderline_hits:
            return "non_chimeric", "Borderline match => rescued"

        # multiple distinct taxonomy => chimeric
        taxa_groups = defaultdict(list)
        for h in db_hits:
            taxa_groups[h["taxonomy"]].append(h)
        if len(taxa_groups) == 1:
            return "borderline", "Single taxonomy, but no borderline or high match"
        else:
            return "chimeric", "Multiple taxonomies, no borderline or high match"

    return "borderline", "Ambiguous"

def write_sequences_to_file(seq_ids, seq_dict, out_fasta):
    """Write seq_ids to out_fasta if non-empty."""
    if not seq_ids:
        logger.info(f"No sequences to write for {out_fasta}.")
        return
    try:
        with open(out_fasta, "w") as fh:
            for sid in sorted(seq_ids):
                if sid in seq_dict:
                    fh.write(f">{sid}\n{seq_dict[sid]}\n")
        logger.info(f"Wrote {len(seq_ids)} sequences to {out_fasta}")
    except Exception as e:
        logger.error(f"Error writing {out_fasta}: {e}")

def write_sequence_details(details, analysis_dir, base_fasta):
    """
    Append a CSV with classification details for each hit.
    """
    if not details:
        logger.debug("No sequence details to record.")
        return

    base_name = os.path.splitext(os.path.basename(base_fasta))[0]
    csv_file = os.path.join(analysis_dir, f"{base_name}_sequence_details.csv")
    need_header = not os.path.isfile(csv_file)

    try:
        with open(csv_file, "a", newline='') as cf:
            writer = csv.writer(cf)
            if need_header:
                writer.writerow([
                    "Sequence ID",
                    "Query Coverage (%)",
                    "Identity (%)",
                    "Classification",
                    "Hit #",
                    "Hit Origin",
                    "Taxonomy"
                ])
            for row in details:
                writer.writerow(row)
        logger.debug(f"Wrote {len(details)} details => {csv_file}")
    except Exception as e:
        logger.error(f"Cannot write details to CSV: {e}")

def parse_blast_results(args):
    """
    Parse a single BLAST XML file, classify each sequence, and write 
    categories to FASTA in temp or analysis directories.
    """
    (xml_path,
     input_chimeras_dir,
     temp_dir,
     analysis_dir,
     high_identity_threshold,
     high_coverage_threshold,
     borderline_identity_threshold,
     borderline_coverage_threshold) = args

    logger.info(f"Starting processing for {xml_path}")
    base_xml = os.path.basename(xml_path)
    # e.g. "ERR6454463.chimeras_blast_results.xml"
    base_name = base_xml.replace("_blast_results.xml", "")

    # We look for chimeras FASTA file with various extensions
    fasta_extensions = [".fasta", ".fa", ".fas", ".fna"]
    fasta_path = None
    
    # Try to find the chimeras file with different extensions
    for ext in fasta_extensions:
        potential_path = os.path.join(input_chimeras_dir, f"{base_name}{ext}")
        if os.path.isfile(potential_path):
            fasta_path = potential_path
            break
    
    # If not found, try with .chimeras suffix
    if not fasta_path:
        for ext in fasta_extensions:
            potential_path = os.path.join(input_chimeras_dir, f"{base_name}.chimeras{ext}")
            if os.path.isfile(potential_path):
                fasta_path = potential_path
                break
    
    if not fasta_path:
        logger.warning(
            f"FASTA file for {xml_path} not found with base name '{base_name}' and extensions {fasta_extensions}. Skipping."
        )
        return set(), set(), set(), set()

    try:
        seqs = load_fasta_sequences(fasta_path)
    except Exception as e:
        logger.warning(f"Cannot load FASTA {fasta_path}: {e}")
        return set(), set(), set(), set()

    non_chimeric = set()
    chimeric = set()
    borderline = set()
    multiple = set()

    details = []

    try:
        with open(xml_path) as handle:
            records = NCBIXML.parse(handle)
            for record in records:
                qid = extract_query_id(record.query)
                if qid not in seqs:
                    continue
                hits_info = analyze_blast_hits(record, qid)
                is_mult = any(h["multiple_hsps"] for h in hits_info)

                if is_mult:
                    multiple.add(qid)
                else:
                    classification, reason = classify_sequence(
                        hits_info,
                        high_identity_threshold,
                        high_coverage_threshold,
                        borderline_identity_threshold,
                        borderline_coverage_threshold
                    )

                    if classification == "non_chimeric":
                        non_chimeric.add(qid)
                    elif classification == "chimeric":
                        chimeric.add(qid)
                    else:  # borderline
                        borderline.add(qid)

                    # Record details for CSV
                    for i, hit in enumerate(hits_info, 1):
                        origin = "Self-sample" if hit["is_same_sample"] else "Database"
                        details.append([
                            qid,
                            f"{hit['coverage']:.2f}",
                            f"{hit['identity']:.2f}",
                            classification,
                            i,
                            origin,
                            hit["taxonomy"]
                        ])

    except Exception as e:
        logger.error(f"Error parsing {xml_path}: {e}")
        return set(), set(), set(), set()

    # Write categorized FASTA files
    # Final naming: <base_prefix>_non_chimeric.fasta, etc.
    # e.g. "ERR6454463_non_chimeric.fasta"
    trimmed_base_name = base_name.replace(".chimeras", "")  # So if base_name=ERR6454462.chimeras => ERR6454462

    if multiple:
        ma_file = os.path.join(analysis_dir, f"{trimmed_base_name}_multiple_alignments.fasta")
        write_sequences_to_file(multiple, seqs, ma_file)

    if non_chimeric:
        nc_file = os.path.join(temp_dir, f"{trimmed_base_name}_non_chimeric.fasta")
        write_sequences_to_file(non_chimeric, seqs, nc_file)

    if borderline:
        bd_file = os.path.join(temp_dir, f"{trimmed_base_name}_borderline.fasta")
        write_sequences_to_file(borderline, seqs, bd_file)

    if chimeric:
        ch_file = os.path.join(analysis_dir, f"{trimmed_base_name}_chimeric.fasta")
        write_sequences_to_file(chimeric, seqs, ch_file)

    if details:
        write_sequence_details(details, analysis_dir, fasta_path)

    logger.info(f"Completed processing for {xml_path}")
    return non_chimeric, chimeric, borderline, multiple

def generate_report(
    all_non_chimeric,
    all_chimeric,
    all_borderline,
    all_multiple,
    file_results,
    output_dir
):
    """Write a final text summary of all classifications."""
    logger.info("=== Generating Chimera Recovery Report ===")

    counts = {
        "Non-Chimeric Sequences": len(all_non_chimeric),
        "Chimeric Sequences": len(all_chimeric),
        "Borderline Sequences": len(all_borderline),
        "Multiple Alignment Sequences": len(all_multiple)
    }
    total = sum(counts.values())
    rep_path = os.path.join(output_dir, "chimera_recovery_report.txt")

    try:
        with open(rep_path, "w") as rf:
            rf.write("Chimera Recovery Report\n")
            rf.write("======================\n\n")
            rf.write("Overall Summary:\n")
            rf.write("----------------\n")
            for cat, val in counts.items():
                rf.write(f"{cat}: {val} ({val/total*100:.1f}%)\n")
            rf.write(f"\nTotal Sequences Processed: {total}\n\n")

            rf.write("Detailed Results by File:\n")
            rf.write("-------------------------\n")
            for xmlf, dct in sorted(file_results.items()):
                rf.write(f"\n{xmlf}:\n")
                for cat, val in dct.items():
                    rf.write(f"  {cat}: {val}\n")
        logger.info(f"Report generated at {rep_path}\n")
    except Exception as e:
        logger.error(f"Error generating report: {e}")

def process_blast_xml_results(
    input_chimeras_dir,
    output_dir,
    high_identity_threshold,
    high_coverage_threshold,
    borderline_identity_threshold,
    borderline_coverage_threshold
):
    """
    Collect all _blast_results.xml in output_dir/xml, parse them,
    and produce final FASTA + text summary.
    """
    logger.info("=== Step 4: Chimera Detection & Recovery ===")

    temp_dir = os.path.join(output_dir, "temp")
    analysis_dir = os.path.join(output_dir, "analysis")
    xml_dir = os.path.join(output_dir, "xml")

    for d in [temp_dir, analysis_dir]:
        clean_directory(d)

    xml_files = [f for f in os.listdir(xml_dir) if f.endswith("_blast_results.xml")]
    if not xml_files:
        logger.error(f"No *_blast_results.xml files in {xml_dir}. Nothing to parse.")
        sys.exit(1)

    file_results = defaultdict(lambda: defaultdict(int))
    all_non_chimeric = set()
    all_chimeric = set()
    all_borderline = set()
    all_multiple = set()

    start_time = time.time()
    log_system_usage()

    args_list = []
    for xmlf in xml_files:
        xml_path = os.path.join(xml_dir, xmlf)
        args_list.append((
            xml_path,
            input_chimeras_dir,
            temp_dir,
            analysis_dir,
            high_identity_threshold,
            high_coverage_threshold,
            borderline_identity_threshold,
            borderline_coverage_threshold
        ))

    n_procs = max(1, multiprocessing.cpu_count() - 1)
    logger.info(f"Using {n_procs} processes for classification.\n")

    with multiprocessing.Pool(processes=n_procs) as pool:
        results = pool.map(parse_blast_results, args_list, chunksize=max(1, len(args_list)//n_procs))

    for i, (nc, ch, bd, ml) in enumerate(results):
        xml_file = args_list[i][0]
        file_results[xml_file] = {
            "Non-Chimeric Sequences": len(nc),
            "Chimeric Sequences": len(ch),
            "Borderline Sequences": len(bd),
            "Multiple Alignment Sequences": len(ml)
        }
        all_non_chimeric.update(nc)
        all_chimeric.update(ch)
        all_borderline.update(bd)
        all_multiple.update(ml)

    # Copy rescued (non-chimeric + borderline) from temp -> output_dir
    for fname in os.listdir(temp_dir):
        if fname.endswith("_non_chimeric.fasta") or fname.endswith("_borderline.fasta"):
            src = os.path.join(temp_dir, fname)
            dst = os.path.join(output_dir, fname)
            shutil.copy(src, dst)
            logger.info(f"Copied {fname} => {dst}")

    # Generate final report
    generate_report(
        all_non_chimeric,
        all_chimeric,
        all_borderline,
        all_multiple,
        file_results,
        output_dir
    )

    log_system_usage()

    # Remove temp_dir
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            logger.info(f"Removed {temp_dir}.")
        except Exception as e:
            logger.error(f"Failed to remove {temp_dir}: {e}")

    total_time = time.time() - start_time
    logger.info(f"Total time for chimera detection: {total_time:.2f}s\n")


###############################################################################
#                                 main()                                      #
###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="False positive chimera detection & recovery for eDNA/Metabarcoding",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--input_chimeras_dir", default="./",
                        help="Directory containing .chimeras.fasta files.")
    parser.add_argument("--self_fasta_dir", default="./",
                        help="Directory with FASTA files for building self-databases.")
    parser.add_argument("--reference_db",
                        help="Path to a reference DB prefix or reference FASTA (script will create DB if needed).")
    parser.add_argument("--output_dir", default="./rescued_reads",
                        help="Directory where all results are written.")
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of CPU threads for BLAST.")
    parser.add_argument("--high_coverage_threshold", type=float, default=99.0,
                        help="High coverage threshold for non-chimeric classification.")
    parser.add_argument("--high_identity_threshold", type=float, default=99.0,
                        help="High identity threshold for non-chimeric classification.")
    parser.add_argument("--borderline_coverage_threshold", type=float, default=89.0,
                        help="Coverage threshold for borderline => non-chimeric.")
    parser.add_argument("--borderline_identity_threshold", type=float, default=80.0,
                        help="Identity threshold for borderline => non-chimeric.")

    args = parser.parse_args()

    # 1. Create Self-Databases
    db_path = os.path.join(args.output_dir, "databases")
    create_self_databases(args.self_fasta_dir, db_path)

    # 2. Reference DB setup
    ref_db_prefix = handle_reference_db(args.reference_db, args.output_dir)
    if ref_db_prefix:
        logger.info(f"Reference database ready at: {ref_db_prefix}")
    else:
        logger.info("No reference database will be used.")

    # 3. Run BLAST
    run_blast_analysis(
        input_chimeras_dir=args.input_chimeras_dir,
        self_db_dir=db_path,
        reference_db_prefix=ref_db_prefix,
        output_dir=args.output_dir,
        threads=args.threads
    )

    # 4. Chimera Detection
    process_blast_xml_results(
        input_chimeras_dir=args.input_chimeras_dir,
        output_dir=args.output_dir,
        high_identity_threshold=args.high_identity_threshold,
        high_coverage_threshold=args.high_coverage_threshold,
        borderline_identity_threshold=args.borderline_identity_threshold,
        borderline_coverage_threshold=args.borderline_coverage_threshold
    )

    logger.info("Pipeline complete. Check your --output_dir for results.\n")


if __name__ == "__main__":
    main()

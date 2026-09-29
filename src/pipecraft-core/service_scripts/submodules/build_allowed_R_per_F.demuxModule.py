#!/usr/bin/python

# Build per-forward-index reverse-index FASTA files for two-step dual-index demux.
# Each unique F index gets a FASTA of only the R indexes (and sample names)
# that are paired with that F in the indexes file.

from Bio import SeqIO
from collections import defaultdict
from sys import argv
import os
import sys

print("python module running: build allowed R indexes per unique F index")

# usage: script.py indexes_file unique_fwd_fasta outdir search_window
if len(argv) != 5:
    sys.exit("usage: build_allowed_R_per_F.demuxModule.py indexes_file unique_fwd_fasta outdir search_window")

script, indexes_file, unique_fwd, outdir, search_window = argv

seq_to_fname = {}
for rec in SeqIO.parse(unique_fwd, "fasta"):
    seq_to_fname[str(rec.seq).upper()] = rec.id

if not seq_to_fname:
    sys.exit("ERROR]: no unique forward indexes found in %s" % unique_fwd)

allowed = defaultdict(list)
for rec in SeqIO.parse(indexes_file, "fasta"):
    parts = str(rec.seq).upper().split("...")
    if len(parts) != 2:
        sys.exit("ERROR]: sample '%s' is not a dual-index record (expected FWD...REV)" % rec.id)
    fseq, rseq = parts[0], parts[1]
    fname = seq_to_fname.get(fseq)
    if fname is None:
        sys.exit("ERROR]: forward index of sample '%s' was not found in unique F file" % rec.id)
    allowed[fname].append((rec.id, rseq))

os.makedirs(outdir, exist_ok=True)
written = 0
for fname, pairs in allowed.items():
    path = os.path.join(outdir, "R_for_%s.fasta" % fname)
    with open(path, "w") as out:
        for sample, rseq in pairs:
            out.write(">%s\nXN{%s}%s\n" % (sample, search_window, rseq))
    written += 1

print("python module finished: wrote %s per-F reverse-index file(s)" % written)

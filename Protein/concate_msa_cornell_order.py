#!/usr/bin/env python3

import sys
import glob
import os
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -----------------------------
# Input directory
# -----------------------------
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <msa_directory>")
    sys.exit(1)

INPUT = sys.argv[1]

# -----------------------------
# Desired gene order
# -----------------------------
gene_order = ["atpD", "recA", "trpB", "gyrB", "rpoB"]

# -----------------------------
# Find files in strict gene order
# -----------------------------
fileList = []

for gene in gene_order:
    matches = glob.glob(os.path.join(INPUT, f"*{gene}*"))
    if not matches:
        print(f"ERROR: No alignment file found for gene '{gene}'")
        sys.exit(1)
    if len(matches) > 1:
        print(f"WARNING: Multiple files found for gene '{gene}', using: {matches[0]}")
    fileList.append(matches[0])

print("Gene order used for concatenation:")
for f in fileList:
    print("  ", os.path.basename(f))

# -----------------------------
# Concatenate by sequence ID
# -----------------------------
concat_seqs = OrderedDict()
gene_lengths = OrderedDict()

for idx, msa in enumerate(fileList):
    gene_name = gene_order[idx]
    current_seqs = {}

    for record in SeqIO.parse(msa, "fasta"):
        current_seqs[record.id] = str(record.seq)

    if idx == 0:
        # Initialize with first gene
        for sid, seq in current_seqs.items():
            concat_seqs[sid] = seq
        gene_lengths[gene_name] = len(next(iter(current_seqs.values())))
    else:
        # Check that IDs match exactly
        if set(current_seqs.keys()) != set(concat_seqs.keys()):
            missing = set(concat_seqs.keys()) ^ set(current_seqs.keys())
            print(f"ERROR: Sequence ID mismatch in {msa}")
            print("Problematic IDs:", ", ".join(missing))
            sys.exit(1)

        for sid in concat_seqs:
            concat_seqs[sid] += current_seqs[sid]

        gene_lengths[gene_name] = len(next(iter(current_seqs.values())))

    print(f"Added {gene_name} ({gene_lengths[gene_name]} bp)")

# -----------------------------
# Write concatenated FASTA
# -----------------------------
with open("merged.fasta", "w") as out:
    for sid, seq in concat_seqs.items():
        record = SeqRecord(Seq(seq), id=sid, description="")
        SeqIO.write(record, out, "fasta")

print("\nConcatenation complete!")
print(f"Total sequences: {len(concat_seqs)}")
print(f"Total alignment length: {len(next(iter(concat_seqs.values())))} bp")
print("Output file: merged.fasta")

# -----------------------------
# OPTIONAL: write partition file
# -----------------------------
with open("partitions.txt", "w") as pf:
    start = 1
    for gene, length in gene_lengths.items():
        end = start + length - 1
        pf.write(f"DNA, {gene} = {start}-{end}\n")
        start = end + 1

print("Partition file written: partitions.txt")

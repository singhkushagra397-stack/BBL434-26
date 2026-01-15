#!/usr/bin/env python3

import sys
from Bio import SeqIO
from collections import defaultdict

# -------------------------
# PARAMETERS
# -------------------------
K = 8       # k-mer length
T = 3       # minimum occurrences
L = 1000    # window length

# -------------------------
# INPUT CHECK
# -------------------------
if len(sys.argv) < 2:
    print("Usage: python script2.py genomic.fa")
    sys.exit(1)

fasta_file = sys.argv[1]

# -------------------------
# READ GENOME
# -------------------------
genome = ""
for record in SeqIO.parse(fasta_file, "fasta"):
    genome += str(record.seq).upper()

n = len(genome)
print(f"Genome length: {n}")

# -------------------------
# FUNCTION TO FIND CLUMPS
# -------------------------
def find_clumps(genome, k, t, L):
    clumps = set()
    counts = defaultdict(int)

    # initialize first window
    first_window = genome[0:L]
    for i in range(len(first_window) - k + 1):
        kmer = first_window[i:i+k]
        counts[kmer] += 1
        if counts[kmer] >= t:
            clumps.add(kmer)

    # slide window across genome
    for start in range(1, n - L + 1):
        end = start + L

        # remove k-mer leaving window
        old_kmer = genome[start - 1 : start - 1 + k]
        if old_kmer in counts:
            counts[old_kmer] -= 1

        # add k-mer entering window
        new_kmer = genome[end - k : end]
        counts[new_kmer] += 1
        if counts[new_kmer] >= t:
            clumps.add(new_kmer)

    return clumps


# -------------------------
# RUN CLUMP FINDER
# -------------------------
clump_kmers = find_clumps(genome, K, T, L)

# -------------------------
# OUTPUT
# -------------------------
print(f"\nDetected (L={L}, k={K}, t={T}) clumps:")
for mer in sorted(clump_kmers):
    print(mer)

print(f"\nTotal clump-forming k-mers: {len(clump_kmers)}")

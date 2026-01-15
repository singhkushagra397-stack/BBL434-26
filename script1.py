#!/usr/bin/env python3

import sys
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

# -------------------------
# PARAMETERS
# -------------------------
K = 8
WINDOW = 5000
STEP = 500

if len(sys.argv) < 2:
    print("Usage: python script1.py genomic.fa")
    sys.exit(1)

fasta_file = sys.argv[1]

# -------------------------
# READ GENOME
# -------------------------
genome = ""
for record in SeqIO.parse(fasta_file, "fasta"):
    genome += str(record.seq).upper()

genome_len = len(genome)
print(f"Genome length: {genome_len}")

# -------------------------
# GLOBAL K-MER COUNTS
# -------------------------
def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

global_counts = Counter(get_kmers(genome, K))

# Avoid divide-by-zero issues
for mer in list(global_counts.keys()):
    if global_counts[mer] == 0:
        global_counts[mer] = 1

# -------------------------
# SLIDING WINDOW
# -------------------------
window_positions = []
enrichment_scores = []

for start in range(0, genome_len - WINDOW + 1, STEP):
    end = start + WINDOW
    window_seq = genome[start:end]

    window_counts = Counter(get_kmers(window_seq, K))

    enrich_vals = []
    for mer, wc in window_counts.items():
        gc = global_counts.get(mer, 1)
        enrich_vals.append(wc / gc)

    avg_enrichment = sum(enrich_vals) / len(enrich_vals)

    window_positions.append(start)
    enrichment_scores.append(avg_enrichment)

print("Enrichment calculation complete.")

# -------------------------
# PLOT AND SAVE
# -------------------------
plt.figure(figsize=(12,4))
plt.plot(window_positions, enrichment_scores)
plt.xlabel("Genome Position")
plt.ylabel(f"Average Enrichment (k={K})")
plt.title("K-mer Enrichment Across Genome")
plt.tight_layout()

plt.savefig("script1_plot.png", dpi=300)
print("Plot saved as script1_plot.png")

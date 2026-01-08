import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python kmer3.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]
k = 3  # fixed k-mer length

for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    kmers = set()

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmers.add(kmer)

    print(f"{record.id}\tUnique_{k}-mers: {len(kmers)}")

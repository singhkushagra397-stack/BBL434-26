import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python script1.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"{record.id}\t{len(record.seq)}")


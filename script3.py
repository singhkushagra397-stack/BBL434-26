import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python count_records.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

count = 0
for _ in SeqIO.parse(fasta_file, "fasta"):
    count += 1

print(f"Number of records: {count}")

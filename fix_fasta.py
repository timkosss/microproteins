import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python check_fasta.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

with open(fasta_file) as handle:
    for i, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
        seq = str(record.seq)
        if not seq:
            print(f"{i}: {record.id} is EMPTY")
        elif any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq):
            print(f"{i}: {record.id} has non-standard characters {set(seq) - set('ACDEFGHIKLMNPQRSTVWY')}")
        elif len(seq) > 2000:
            print(f"{i}: {record.id} is very long ({len(seq)})")
    print("Done")
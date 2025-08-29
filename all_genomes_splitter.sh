#!/bin/bash
input="/wistar/auslander/Timothy/microproteins_timka/decr_thresh_bad_proteins_0.6_merged.fa"
outdir="/wistar/auslander/Timothy/microproteins_timka/bad_protein_chunks"
mkdir -p "$outdir"

# Count sequences and calculate per-chunk size
total=$(grep -c "^>" "$input")
per_chunk=$(( (total + 9) / 10 ))

awk -v per_chunk=$per_chunk -v outdir="$outdir" -v base="db2_chunk" '
  BEGIN {chunk=1; seq=0; outfile=sprintf("%s/%s_%02d.fa", outdir, base, chunk)}
  /^>/ {
    if (seq >= per_chunk) {
      chunk++;
      seq=0;
      outfile=sprintf("%s/%s_%02d.fa", outdir, base, chunk);
    }
    seq++;
  }
  {print > outfile}
' "$input"

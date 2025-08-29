#!/bin/bash
#SBATCH --job-name=merge_cdhit2d
#SBATCH --output=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/merge_cdhit2d.out
#SBATCH --error=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/merge_cdhit2d.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

OUTDIR="/wistar/auslander/Timothy/microproteins_timka/bad_protein_chunk_outputs"
FINAL="/wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_merged.fa"
FINAL_CLSTR="/wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_merged.fa.clstr"

echo "Merging outputs from $OUTDIR"

# Merge FASTA outputs
cat $OUTDIR/bad_protein_chunk_*_vs_prodigal.fa > "$FINAL"

# Merge cluster outputs (optional)
cat $OUTDIR/bad_protein_chunk_*_vs__vs_prodigal.clstr > "$FINAL_CLSTR"

echo "Final merged FASTA: $FINAL"
echo "Final merged CLSTR: $FINAL_CLSTR"

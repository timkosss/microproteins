#!/bin/bash
#SBATCH --job-name=merge_cdhit2d
#SBATCH --output=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/merge_cdhit2d.out
#SBATCH --error=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/merge_cdhit2d.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

OUTDIR="/wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/negative_protein_chunks/"
FINAL="/wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/super_negative_proteins.fa"


echo "Merging outputs from $OUTDIR"

# Merge FASTA outputs
cat $OUTDIR/* > "$FINAL"

# Merge cluster outputs (optional)

echo "Final merged FASTA: $FINAL"

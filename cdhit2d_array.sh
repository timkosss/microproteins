#!/bin/bash
#SBATCH --job-name=cdhit2d_array
#SBATCH --output=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/cdhit2d_%A_%a.out
#SBATCH --error=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/cdhit2d_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-10
source /applications/miniforge3/23.11.4/etc/profile.d/conda.sh
conda activate mycdhit
# Paths
DB1="/wistar/auslander/Timothy/microproteins_timka/complete_prodigal.fa"
CHUNK_DIR="/wistar/auslander/Timothy/microproteins_timka/bad_protein_chunks"
OUTDIR="/wistar/auslander/Timothy/microproteins_timka/bad_protein_chunk_outputs"
mkdir -p "$OUTDIR"

# Map array index -> chunk file
CHUNK=$(printf "%s/db2_chunk_%02d.fa" "$CHUNK_DIR" "$SLURM_ARRAY_TASK_ID")

# Output file
OUTFILE="$OUTDIR/bad_protein_chunk_${SLURM_ARRAY_TASK_ID}_vs_prodigal.fa"

echo "Processing $CHUNK"
cd-hit-2d \
  -i "$DB1" \
  -i2 "$CHUNK" \
  -o "$OUTFILE" \
  -c 0.6 -n 4 -M 50000 -T 8

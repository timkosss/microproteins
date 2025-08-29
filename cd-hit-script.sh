#!/bin/bash
#SBATCH --job-name=cdhit2d_job_0.6_64cpu
#SBATCH --output=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/cdhit2d_%j.out
#SBATCH --error=/wistar/auslander/Timothy/microproteins_timka/OUTFiles/cdhit2d_%j.err
#SBATCH --time=2-00:00:00          # max runtime
#SBATCH --mem=100G                # memory requirement
#SBATCH --cpus-per-task=44        # threads

source /applications/miniforge3/23.11.4/etc/profile.d/conda.sh
conda activate mycdhit

# Input files
DB1="/wistar/auslander/Timothy/microproteins_timka/microprotein_fixed.fa"
DB2="/wistar/auslander/Timothy/microproteins_timka/all_genomes_protein.fasta"
DB3="/wistar/auslander/Timothy/microproteins_timka/complete_prodigal.fa"
OUT1="/wistar/auslander/Timothy/microproteins_timka/decr_thresh_bad_proteins_0.6.fa"
OUT2="/wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6.fa"

# Run CD-HIT-2D
#cd-hit-2d -i $DB1 -i2 $DB2 -o $OUT1 -c 0.6 -n 4 -M 50000 -T 8

#cd-hit-2d -i $DB3 -i2 $OUT1 -o $OUT2 -c 0.6 -n 4 -M 50000 -T 8

# cd-hit \
# -i /wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_merged.fa \
# -o /wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_clustered.fa \
# -c 0.6 \
# -n 4 \
# -M 50000 \
# -T 16

cd-hit \
-i /wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_merged.fa \
-o /wistar/auslander/Timothy/microproteins_timka/decr_thresh_super_bad_proteins_0.6_clustered2.fa \
-c 0.6 \
-n 4 \
-M 50000 \
-T 64


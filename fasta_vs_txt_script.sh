#!/bin/bash

#SBATCH -J compare_fasta_vs_txt
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/outfiles/output_%A_%a.out
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/outfiles/output_%A_%a.err
#SBATCH --mem=1G
#SBATCH --array=1-5667
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00

echo "Running On:"
srun hostname
srun uname -a

FILE=$(awk "NR==$SLURM_ARRAY_TASK_ID" /wistar/auslander/microproteins/splitted_fa/sample.txt)

echo "Starting Python Script on $MAIN_FILE"
python -u /wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/is_in_prodigal.py "$FILE"

#################################
#           end   script        #
#################################
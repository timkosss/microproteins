#!/bin/bash

#SBATCH -J create_fasta_with_loc
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mem=1G
#SBATCH --array=1-5667
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00

echo "Running On:"
srun hostname
srun uname -a

FILE=$(awk "NR==$SLURM_ARRAY_TASK_ID" /wistar/auslander/microproteins/splitted_fa/sample.txt)

MAIN_FILE="/wistar/auslander/microproteins/splitted_fa/${FILE}"

echo "Starting Python Script on $MAIN_FILE"
python -u /wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/valid_protein_loc.py "$MAIN_FILE"

#################################
#           end   script        #
#################################
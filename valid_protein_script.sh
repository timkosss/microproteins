#!/bin/bash

#SBATCH -J make_protein_fasta
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/OUTFiles/output_%A.txt
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/OUTFiles/output_%A.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00

echo "Running On:"
srun hostname
srun uname -a

echo "Starting Python Script..."
python -u /wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/valid_protein_loc.py
#################################
#           end   script        #
#################################

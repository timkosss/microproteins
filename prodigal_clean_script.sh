#!/bin/bash

#SBATCH -J prodigal_cleaning
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/prodigal_OUTFIles/output2.txt
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/prodigal_OUTFIles/output2.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00

echo "Running On:"
srun hostname
srun uname -a

echo "Starting Python Script..."
python -u /wistar/auslander/Timothy/microproteins_timka/prodigal_to_fasta.py 

#################################
#           end   script        #
#################################

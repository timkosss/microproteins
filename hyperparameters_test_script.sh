#!/bin/bash

#SBATCH -J hyperparameters_test
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/OUTFiles/hyperparameters_output_%A_%a.txt
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/OUTFiles/hyperparameters_output_%A_%a.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00
#SBATCH --array=1-7

echo "Running On:"
srun hostname
srun uname -a

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    echo "AUROC,LOSS,SEQ_LENGTH,EPOCH,LRS,DB,TRAINING,DROPOUT" > hyperparameter_results_negative_num_training.csv 
fi

lrs=(0.0001 0.0002 0.0005 0.001 0.002 0.005)

seq_lens=(10 15 20 25)
dropouts=(0.4 0.5 0.6)

epochs=(15 20 25)
neg_dbs=("contained_within_prodigal/super_negative_proteins.fa" "intermediate_files/decr_thresh_super_bad_proteins_0.6_clustered2.fa")
num_training=(1000 5000 10000 20000 25000 30000 32500)

training=${num_training[($SLURM_ARRAY_TASK_ID-1)]}
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"
python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "15" "0.001" "contained_within_prodigal/super_negative_proteins.fa" "$training" "20" "0.45"



# lr=${lrs[($SLURM_ARRAY_TASK_ID-1)]}
# for seq_len in ${seq_lens[@]}; do
#     for dropout in ${dropouts[@]}; do
#         for epoch in ${epochs[@]}; do
#             for db in ${neg_dbs[@]}; do
#                 for training in ${num_training[@]}; do
#                     echo "Starting Python Script..."
#                     python -u /wistar/auslander/Timothy/microproteins_timka/microprot_timka.py "$epoch" "$lr" "$db" "$training" "$seq_len" "$dropout"
#                 done 
#             done 
#         done
#     done
# done 

#################################
#           end   script        #
#################################

#!/bin/bash

#SBATCH --job-name=pinn_vlasov
#SBATCH --nodes=1                    # we request one node
#SBATCH --ntasks-per-node=1          # with one task per node (= number of GPUs here)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
##SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --partition=gpu_p2l
#SBATCH -o pinn_voice.out
#SBATCH -e pinn_voice.err
#SBATCH --qos=qos_gpu-dev
#SBATCH -A wjo@gpu
set +v

arg=$#

if [ $# -lt 2 ]; then
    echo " Syntaxe:  jean-zay-gpu-pinn.slurm <YAML_input_PINN> <VOICEXX_directory>"
    exit
fi

##----- Configuration de l'environnement -----
module load tensorflow-gpu

python PINNmain.py -i $1 -p $2

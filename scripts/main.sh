#!/bin/bash

#SBATCH -p evlab
#SBATCH --ntasks=1 
#SBATCH --time=3:00:00
#SBATCH --mem=4G
#SBATCH --output=./slurm_log/main%j.out

source ~/.local/bin/activate-conda
conda activate r_glotto

echo 'executing main.sh'
date

Rscript main.R


echo 'finished!'
date
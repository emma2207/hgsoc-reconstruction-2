#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mem=128M
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=split_barcodes
#SBATCH --output=split_barcodes_%J.log
#SBATCH --error=split_barcodes_%J.err

module load miniforge
conda activate python_vireosnp_install

# Script location is /scratch/alpine/$USER/hgsoc/demultiplexing
python split_barcodes.py

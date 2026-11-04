#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=16
#SBATCH --account=amc-general
#SBATCH --job-name=cellranger_multi
#SBATCH --output=01_cellranger_multi_%A_%a.log
#SBATCH --error=01_cellranger_multi_%A_%a.err

module use --append /projects/$USER/lmod-files
module load cellranger/9.0.1

export pool=$SLURM_ARRAY_TASK_ID

# Location is /scratch/alpine/$USER/hgsoc/demultiplexing/cellranger

cellranger multi --id=pool$pool \
	--csv=cellranger_input/cellranger_multi_input_pool$pool.csv

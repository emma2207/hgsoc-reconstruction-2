#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=htseq-count
#SBATCH --output=counts_%A_%a.log
#SBATCH --error=htseq-count_%A_%a.err

module load anaconda
conda activate htseq_install

export pool=$SLURM_ARRAY_TASK_ID # 1-10
sample=$1                        # e.g. 2018

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
original_location=$(pwd)
index_location="/projects/$USER/hgsoc/refdata-gex-GRCh38-2024-A"
sample_location="$original_location/subset-bam/pool$pool/HGSOC-$sample"

htseq-count -r "name" \
	-f "bam" \
	$sample_location/sample.bam \
	$index_location/genes/genes.gtf

#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=cellsnp-lite
#SBATCH --output=03_genotype_pooled_%A_%a.log
#SBATCH --error=03_genotype_pooled_%A_%a.err

module load anaconda
conda activate cellsnp-lite_install

export pool=$SLURM_ARRAY_TASK_ID

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
mkdir -p results/cellSNP_no_ref/pool${pool}

original_location=$(pwd)
pooled_bam_location="$original_location/bam/pool${pool}"
barcodes_zip=$(ls /pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC/pooled_sc/Pool${pool}-GEX-CaseyGreene-*/outs/raw_feature_bc_matrix/barcodes.tsv.gz)
output_location="$original_location/results/cellSNP_no_ref/pool${pool}"

cellsnp-lite \
	-s $pooled_bam_location/pooled.bam \
	-b "${barcodes_zip}" \
	-O $output_location \
	-p 20 \
	--minMAF=0.1 \
	--minCOUNT=100 \
	--gzip

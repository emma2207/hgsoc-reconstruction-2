#!/bin/bash

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=split-bam
#SBATCH --output=split-bam_%A_%a.log
#SBATCH --error=split-bam_%A_%a.err

module load miniforge
conda activate subset-bam-install

export pool=$SLURM_ARRAY_TASK_ID

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
original_location=`pwd`
sample_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/aligned_reads/hgsoc-new/pooled_single_cell/pool${pool}"

find "htseq-count/pool${pool}/" -type d -name "HGSOC-[0-9][0-9][0-9][0-9]" -print0 | while IFS= read -r -d $'\0' subdir; do
# find "htseq-count/pool${pool}/" -type d -name "donor*" -print0 | while IFS= read -r -d $'\0' subdir; do
    echo "Processing subdirectory: ${subdir}"
    barcodes_location=${subdir}
    sample=$(basename "${subdir}")
    echo "Processing sample: ${sample}"
    mkdir -p subset-bam/pool${pool}/${sample}

subset-bam -b ${sample_location}/pooled.bam \
    -c ${barcodes_location}/barcodes_per_sample.csv \
    -o "subset-bam/pool${pool}/${sample}/sample.bam" \
    --cores 1
done

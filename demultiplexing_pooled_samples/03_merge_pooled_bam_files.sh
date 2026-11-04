#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=merge_pooled_bam
#SBATCH --output=merge_pooled_bam_%A_%a.log
#SBATCH --error=merge_pooled_bam_%A_%a.err

module load samtools/1.16.1

export pool=$SLURM_ARRAY_TASK_ID

# Script location is /scratch/alpine/$USER/hgsoc/demultiplexing
original_location=$(pwd)
sample_alignment_location="$original_location/cellranger/pool$pool/outs/per_sample_outs"
unassigned_location="$original_location/cellranger/pool$pool/outs/multi/count"
output_location="$original_location/bam/pool$pool"

samtools merge $output_location/pooled.bam \
	$unassigned_location/unassigned_alignments.bam \
	$sample_alignment_location/sample1/count/sample_alignments.bam \
	$sample_alignment_location/sample2/count/sample_alignments.bam \
	$sample_alignment_location/sample3/count/sample_alignments.bam \
	$sample_alignment_location/sample4/count/sample_alignments.bam

samtools index $output_location/pooled.bam

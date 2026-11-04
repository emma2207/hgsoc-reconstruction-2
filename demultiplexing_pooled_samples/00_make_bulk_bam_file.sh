#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=star
#SBATCH --output=01_make_bam_file_%A_%a.log
#SBATCH --error=01_make_bam_file_%A_%a.err

module load star/2.7.10b
module load samtools/1.16.1

folder=$1  # 230414, 230418, 230509, 230626
export taskID=$SLURM_ARRAY_TASK_ID

# Original location is /scratch/alpine/${USER}/hgsoc/bulk
original_location=`pwd`
index_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A/star"
pl_location="/pl/active/cgreene-sc-hgsoc"

if [[ "${folder}" == "230414" ]]; then
    fastq_location="${pl_location}/ariel_sc_HGSOC/HippenA_230414_A00901_0909_BH3FJ2DRX3"
    config="${original_location}/230414_config.txt"
elif [[ "${folder}" == "230418" ]]; then
    fastq_location="${pl_location}/ariel_sc_HGSOC/HippenA_230418_A00901_0912_BHK72LDMXY"
    config="${original_location}/230418_config.txt"
elif [[ "${folder}" == "230509" ]]; then
    fastq_location="${pl_location}/ariel_sc_HGSOC/HippenA_230509_A00901_0926_AHK3THDMXY"
    config="${original_location}/230509_config.txt"
elif [[ "${folder}" == "230626" ]]; then
    fastq_location="${pl_location}/penn_HGSOC"
    config="${original_location}/230626_config.txt"
else
    echo "Folder ${folder} not found. Check your input."
fi

sample=$(awk -v ArrayTaskID=${taskID} '$1==ArrayTaskID {print $2}' ${config})

mkdir -p ${folder}/${sample}
output_location="${original_location}/${folder}/${sample}"

#STAR only accepts one set of files, so we'll merge across lanes
cat ${fastq_location}/HGSOC-${sample}_*/HGSOC-${sample}_S*_L001_R1_001.fastq.gz \
    ${fastq_location}/HGSOC-${sample}_*/HGSOC-${sample}_S*_L002_R1_001.fastq.gz \
    > ${output_location}/HGSOC-${sample}_R1_merged.fastq.gz
cat ${fastq_location}/HGSOC-${sample}_*/HGSOC-${sample}_S*_L001_R2_001.fastq.gz \
    ${fastq_location}/HGSOC-${sample}_*/HGSOC-${sample}_S*_L002_R2_001.fastq.gz \
    > ${output_location}/HGSOC-${sample}_R2_merged.fastq.gz

STAR \
	--genomeDir ${index_location} \
	--runThreadN 6 \
	--readFilesIn ${output_location}/HGSOC-${sample}_R1_merged.fastq.gz ${output_location}/HGSOC-${sample}_R2_merged.fastq.gz \
	--outFileNamePrefix ${output_location}/STAR/ \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts

samtools index ${output_location}/STAR/Aligned.sortedByCoord.out.bam

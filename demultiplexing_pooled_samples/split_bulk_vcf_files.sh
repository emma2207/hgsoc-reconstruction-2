#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=02_bcftools_%J.log
#SBATCH --error=02_bcftools_%J.err

module load bcftools/1.16

# Script location is /scratch/alpine/$USER/hgsoc/pooled
original_location=`pwd`
data_location="${original_location}/bcftools/all"

input_file_bulk="${data_location}/bcftools_bulk_all_rehead.vcf"
input_file_diss_bulk="${data_location}/bcftools_diss_bulk_all_rehead.vcf"

for sample in `bcftools query -l ${input_file_bulk}` 
do
    short_sample_id="${sample: -4}"

    bcftools view \
        -c1 -Ov -s $sample \
        -o "${data_location}/bcftools_diss_bulk_sample${short_sample_id}.vcf" \
        ${input_file}

    bcftools view \
        -i 'QUAL>20' "${data_location}/bcftools_diss_bulk_sample${short_sample_id}.vcf" \
        > "${data_location}/bcftools_filtered_diss_bulk_sample${short_sample_id}.vcf"

    bcftools view \
        -i 'DP>30' "${data_location}/bcftools_filtered_bulk_sample${short_sample_id}.vcf" \
        > "${data_location}/bcftools_filtered2_bulk_sample${short_sample_id}.vcf"

done

for sample in `bcftools query -l ${input_file_diss_bulk}` 
do
    short_sample_id="${sample: -4}"

    bcftools view \
        -c1 -Ov -s $sample \
        -o "${data_location}/bcftools_diss_bulk_sample${short_sample_id}.vcf" \
        ${input_file}

    bcftools view \
        -i 'QUAL>20' "${data_location}/bcftools_diss_bulk_sample${short_sample_id}.vcf" \
        > "${data_location}/bcftools_filtered_diss_bulk_sample${short_sample_id}.vcf"

    bcftools view \
        -i 'DP>30' "${data_location}/bcftools_filtered_diss_bulk_sample${short_sample_id}.vcf" \
        > "${data_location}/bcftools_filtered2_diss_bulk_sample${short_sample_id}.vcf"

done
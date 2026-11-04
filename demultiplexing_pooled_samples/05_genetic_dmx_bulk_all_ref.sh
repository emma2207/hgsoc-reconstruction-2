#!/bin/bash

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=vireo
#SBATCH --output=04_run_vireo_%A_%a.log
#SBATCH --error=04_run_vireo_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate vireo_install

export pool=$SLURM_ARRAY_TASK_ID
data_type=$1  # "bulk_chunk_ribo" or "diss_bulk_polyA"
samples_to_exclude=$2  # Sample to exclude from the bulk VCF file (comma separated)

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
mkdir -p results/vireo/pool${pool}

original_location=$(pwd)
bulk_ref_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/variant_calls/hgsoc-new"
cellsnp_location="$original_location/results/cellSNP/pool${pool}"
output_location="$original_location/results/vireo/pool${pool}"

# Adjust the number of samples for pools 9 & 10
if [ "$pool" -lt 9 ]; then
	n_samples=4
elif [ "$pool" -lt 10 ]; then
	n_samples=3
else
	n_samples=2
fi

# Select the appropriate bulk VCF file based on data type
if [ "${data_type}" == "diss_bulk_polyA" ]; then 
	bulk_vcf_file=$bulk_ref_location/bulk_diss_polyA_all_samples.vcf.gz
elif [ "${data_type}" == "bulk_chunk_ribo" ]; then 
	bulk_vcf_file=$bulk_ref_location/bulk_chunk_ribo_all_samples.vcf.gz
else
	echo "Invalid data type specified. Use 'bulk_chunk_ribo' or 'diss_bulk_polyA'."
	exit 1
fi

# Optionally exclude low quality samples from the bulk VCF file.
# Each entry in samples_to_exclude is a partial name matched against the VCF sample names.
if [ -z "$samples_to_exclude" ]; then
	echo "No samples to exclude; using the original reference VCF."
	filtered_bulk_vcf_file=$bulk_vcf_file
else
	IFS=',' read -ra samples <<< "$samples_to_exclude"
	full_samples_to_exclude=""
	for sample in "${samples[@]}"; do
		echo "Excluding sample: $sample"
		matches=$(bcftools query -l "$bulk_vcf_file" | grep -E "^HGSOC-(${sample})_" || true)
		if [ -z "$matches" ]; then
			echo "No VCF sample matched '${sample}'. Expected a sample starting with HGSOC-(${sample})_."
			exit 1
		fi
		match_count=$(printf "%s\n" "$matches" | wc -l | tr -d ' ')
		if [ "$match_count" -ne 1 ]; then
			echo "Expected exactly 1 VCF sample match for '${sample}', found ${match_count}:"
			printf "%s\n" "$matches"
			exit 1
		fi
		full_samples_to_exclude+="${matches},"
	done
	full_samples_to_exclude="${full_samples_to_exclude%,}"

	filtered_bulk_vcf_file=${bulk_vcf_file%.vcf.gz}_pool${pool}_rehead.vcf.gz
	bcftools view -s ^${full_samples_to_exclude} "$bulk_vcf_file" -Oz -o "$filtered_bulk_vcf_file"
	bcftools index -t -f "$filtered_bulk_vcf_file"
fi

# Run vireo
vireo \
	-c $cellsnp_location \
	-N $n_samples \
	-o $output_location \
	-d $filtered_bulk_vcf_file \
	--randSeed=12

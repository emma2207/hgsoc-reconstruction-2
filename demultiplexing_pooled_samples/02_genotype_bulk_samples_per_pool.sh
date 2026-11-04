#!/bin/sh

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=5G
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=02_bcftools_%A_%a.log
#SBATCH --error=02_bcftools_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools/1.16

export pool=$SLURM_ARRAY_TASK_ID
data_type=$1  # "diss_bulk_polyA" or "bulk_chunk_ribo"
pool_csv=$2  # CSV file containing samples per pool

# Script location is /scratch/alpine/${USER}/hgsoc/demultiplexing
original_location=`pwd`
index_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A"
bulk_bam_location="/scratch/alpine/${USER}/hgsoc/bulk"
mkdir -p bcftools/pool${pool}
output_location="${original_location}/bcftools/pool${pool}"

# Figure out what samples are in the pool
while IFS=, read -r index col1 col2 col3 col4
do
	if [[ "$index" == "Pool${pool}" ]]; then
		col1=${col1%.0}
		col2=${col2%.0}
		col3=${col3%.0}
		col4=${col4%.0}
		sample_array=($col1 $col2 $col3 $col4)
		echo ${sample_array[*]}
	fi
done < ${pool_csv}

# Check which folder the bulk seq data for each sample is in
# Dissociated bulk samples
samples_230414=("2018" "2126" "2202" "2221" "2240" \
"2278" "2364" "2407" "2416" "2477" "2507" \
"2023" "2129" "2209" "2238" "2249" "2313" \
"2401" "2408" "2460" "2483" "2526")
samples_230418=("2094" "2186" "2216" "2230" "2246" \
"2309" "2423" "2455" "2466" "2514")
# Bulk samples
samples_230509=("2018" "2023" "2094" "2126" "2129" \
	"2186" "2202" "2209" "2240" "2309" "2313" "2401" \
	"2407" "2444" "2483" "2526")
samples_230626=("2216" "2221" "2230" "2238" "2246" \
	"2249" "2268" "2278" "2296" "2364" "2408" "2416" \
	"2423" "2430" "2455" "2460" "2466" "2477" "2507" \
	"2514")

bam_files=()
for sample in "${sample_array[@]}"
do
	if [[ "$data_type" == "diss_bulk_polyA" ]]; then
		if [[ " ${samples_230414[*]} " =~ " ${sample} " ]]; then
			bam_files+=("${bulk_bam_location}/230414/${sample}/STAR/Aligned.sortedByCoord.out.bam")
		elif [[ " ${samples_230418[*]} " =~ " ${sample} " ]]; then
			bam_files+=("${bulk_bam_location}/230418/${sample}/STAR/Aligned.sortedByCoord.out.bam")
		else
			echo "Sample ${sample} not found in dissociated bulk."
		fi
	elif [[ "$data_type" == "bulk_chunk_ribo" ]]; then
		if [[ " ${samples_230509[*]} " =~ " $sample " ]]; then
			bam_files+=("${bulk_bam_location}/230509/${sample}/STAR/Aligned.sortedByCoord.out.bam")
		elif [[ " ${samples_230626[*]} " =~ " $sample " ]]; then
			bam_files+=("${bulk_bam_location}/230626/${sample}/STAR/Aligned.sortedByCoord.out.bam")
		else
			echo "Sample ${sample} not found in bulk."
		fi
	else
		echo "Invalid data type: ${data_type}. Use 'diss_bulk_polyA' or 'bulk_chunk_ribo'."
		exit 1
	fi
done

echo ${bam_files[@]}
nr_bam_files=${#bam_files[@]}
read_depth_cutoff=$((30 * ${nr_bam_files}))
echo "Number of BAM files: ${nr_bam_files}"

# Run bcftools
bcftools mpileup -Ou \
	-f ${index_location}/fasta/genome.fa \
	${bam_files[@]} | \
bcftools call -mv -Ov \
	-o ${output_location}/bcftools_${data_type}_pool${pool}.vcf

bcftools reheader \
	-s bcftools/bcftools_${data_type}_rename.txt \
	${output_location}/bcftools_${data_type}_pool${pool}.vcf > \
	${output_location}/bcftools_${data_type}_pool${pool}_rehead.vcf 

# Filtering with bcftools
# Quality score greater than 20
bcftools view \
	-i "QUAL>20" ${output_location}/bcftools_${data_type}_pool${pool}_rehead.vcf \
	> ${output_location}/bcftools_qual_filter_${data_type}_pool${pool}.vcf
# Read depth greater than 30 x # of samples
bcftools view \
	-i "DP>${read_depth_cutoff}" ${output_location}/bcftools_qual_filter_${data_type}_pool${pool}.vcf \
	> ${output_location}/bcftools_qual_dp_filter_${data_type}_pool${pool}.vcf

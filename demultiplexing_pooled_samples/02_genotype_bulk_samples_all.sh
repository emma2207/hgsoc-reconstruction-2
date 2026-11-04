#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=8
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=02_bcftools_%J.log
#SBATCH --error=02_bcftools_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools/1.16

data_type=$1  # "diss_bulk_polyA" or "bulk_chunk_ribo"

# Script location is /scratch/alpine/${USER}/hgsoc/demultiplexing
original_location=`pwd`
index_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A"
bulk_bam_location="/scratch/alpine/${USER}/hgsoc/bulk"
mkdir -p results/bcftools/all
output_location="${original_location}/results/bcftools/all"

# Check which folder the bulk data for each sample is in
# Dissociated Bulk
samples_230414=("2018" "2126" "2202" "2221" "2240" \
	"2278" "2364" "2407" "2416" "2477" "2507" \
	"2023" "2129" "2209" "2238" "2249" "2313" \
	"2401" "2408" "2460" "2483" "2526")
samples_230418=("2094" "2186" "2216" "2230" "2246" \
	"2309" "2423" "2455" "2466" "2514")
# Bulk
samples_230509=("2018" "2023" "2094" "2126" "2129" \
	"2186" "2202" "2209" "2240" "2309" "2313" "2401" \
	"2407" "2444" "2483" "2526")
samples_230626=("2216" "2221" "2230" "2238" "2246" \
	"2249" "2268" "2278" "2296" "2364" "2408" "2416" \
	"2423" "2430" "2455" "2460" "2466" "2477" "2507" \
	"2514")

bam_files=()
if [ "${data_type}" == "diss_bulk_polyA" ]; then
	# Dissociated Bulk samples
	for sample in ${samples_230414[@]}
	do
		bam_files+=("${bulk_bam_location}/230414/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	done

	for sample in ${samples_230418[@]}
	do
		bam_files+=("${bulk_bam_location}/230418/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	done

elif [ "${data_type}" == "bulk_chunk_ribo" ]; then
	# Bulk samples
	for sample in ${samples_230509[@]}
	do
		bam_files+=("${bulk_bam_location}/230509/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	done

	for sample in ${samples_230626[@]}
	do
		bam_files+=("${bulk_bam_location}/230626/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	done
else
	echo "Invalid data type specified. Use 'diss_bulk_polyA' or 'bulk_chunk_ribo'."
	exit 1
fi

# Genotyping with bcftools
bcftools mpileup -Ou \
	-f ${index_location}/fasta/genome.fa \
	${bam_files[@]} | \
bcftools call -mv -Ov \
	-o ${output_location}/bcftools_${data_type}_all.vcf

bcftools reheader \
	-s bcftools_${data_type}_rename.txt \
	${output_location}/bcftools_${data_type}_all.vcf > \
	${output_location}/bcftools_${data_type}_all_rehead.vcf 

# Filtering with bcftools
# Quality score greater than 20
bcftools view \
	-i 'QUAL>20' ${output_location}/bcftools_${data_type}_all_rehead.vcf \
	> ${output_location}/bcftools_qual_filter_${data_type}_all.vcf
# Read depth greater than 900 (roughly 30 x # of samples)
bcftools view \
	-i 'DP>900' "${output_location}/bcftools_qual_filter_${data_type}_all.vcf" \
	> "${output_location}/bcftools_qual_dp_filter_${data_type}_all.vcf"
# Remove sample 2507 (poor quality sample)
if [ "${data_type}" == "diss_bulk_polyA" ]; then
	bcftools view \
	-s ^"230414/2507" "${output_location}/bcftools_qual_dp_filter_${data_type}_all.vcf" \
	> "${output_location}/bcftools_qual_dp_rm_2507_filter_${data_type}_all.vcf"
fi

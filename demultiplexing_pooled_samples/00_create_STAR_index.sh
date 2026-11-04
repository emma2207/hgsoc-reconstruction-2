#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=star
#SBATCH --output=00_create_STAR_index_%J.log
#SBATCH --error=00_create_STAR_index_%J.err

module load star/2.7.10b

# See README for details about the reference genome.

# Current location is /scratch/alpine/${USER}/hgsoc/scripts
current_location=`pwd`
index_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A"

# Check if STAR index exists
if [ -f "$index_location/star/SAindex" ]; then
	echo "STAR index exists, skipping."
else
	STAR \
		--runMode genomeGenerate \
		--runThreadN 20 \
		--genomeDir $index_location/star \
		--genomeFastaFiles $index_location/fasta/genome.fa \
		--sjdbGTFfile $index_location/genes/genes.gtf
fi

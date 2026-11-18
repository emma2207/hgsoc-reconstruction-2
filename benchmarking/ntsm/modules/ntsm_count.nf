#!/usr/bin/env nextflow

process NTSM_COUNT {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_counts", mode: 'copy'

    input:
    tuple val(modality), val(sample_name), path(r1_fastqs), path(r2_fastqs)

    output:
    path("counts_${modality}_${sample_name}.txt")

    script:
    def sites = file(params.sites).toAbsolutePath()
    """
    set -euo pipefail

    output_datatype="${params.pseudobulk ? 'pseudobulk' : params.datatype}"
    output_location="ntsm_counts/${params.dataset}/\${output_datatype}/${modality}"
    mkdir -p \${output_location}

    cat ${r1_fastqs} > \${output_location}/${sample_name}.R1.fastq.gz
    cat ${r2_fastqs} > \${output_location}/${sample_name}.R2.fastq.gz

    if [ "${params.read_type}" = "single" ]
    then
        ntsmCount -t ${task.cpus} -s "${sites}" \
            \${output_location}/${sample_name}.R1.fastq.gz > counts_${modality}_${sample_name}.txt
    else
        ntsmCount -t ${task.cpus} -s "${sites}" \
            \${output_location}/${sample_name}.R1.fastq.gz \
            \${output_location}/${sample_name}.R2.fastq.gz > counts_${modality}_${sample_name}.txt
    fi
    """
}

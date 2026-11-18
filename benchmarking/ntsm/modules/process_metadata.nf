#!/usr/bin/env nextflow

process PROCESS_METADATA {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/metadata/${params.dataset}/${params.read_type}", mode: 'copy'

    input:
    tuple val(modality), path(metadata), val(fastq_dir)

    output:
    tuple val(modality), path("run_name_mapping_${modality}.csv")

    script:
    def metadataScript = file("${params.projectDir}/process_metadata.py").toAbsolutePath()
    """
    python "${metadataScript}" \
        --metadata \"$metadata\" \
        --fastq-dir \"$fastq_dir\" \
        --read-type ${params.read_type} \
        --output "run_name_mapping_${modality}.csv"
    """
}

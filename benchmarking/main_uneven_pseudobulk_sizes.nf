#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variantion on main_from_genotype.nf that can be used to assess
 * the impact of comparing different size pseudobulks on sample matching performance.
 * This pipeline analyzes VCF files to determine if samples come from the same donor
 * using multiple tools: HYSYS, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

include { validateBenchmarkParameters } from './modules/helpers/modality_utils'
include {
    requirePositiveInteger;
    validateVcfRunPaths
} from './modules/helpers/path_validation'

// Import processes
include { MERGE_AND_FILTER_VCFS } from './modules/1a_merge_vcfs'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'


// Workflow specific parameters
params.ncells_2 = 100000 // Number of cells included in second set of pseudobulk samples (for comparison with params.ncells)
params.outdir = "results/uneven_pseudobulk_sizes/ncells_${params.ncells_2}"


// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)
    if (!params.pseudobulk) {
        throw new IllegalArgumentException("main_uneven_pseudobulk_sizes.nf requires params.pseudobulk = true.")
    }
    requirePositiveInteger(params.ncells_2, 'params.ncells_2')
    if (params.ncells_2 == params.ncells) {
        throw new IllegalArgumentException("params.ncells_2 must differ from params.ncells for an uneven-size comparison.")
    }

    validateVcfRunPaths(
        params,
        ['vireo', 'fingerprints', 'hysys', 'ngscheckmate'],
        [
            set_1: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz",
            set_2: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells_2}/*_2.vcf.gz"
        ]
    )

    // Find all vcfs with the listed datatypes
    modalities = ["pseudobulk"]
    mod_channel = channel.fromList(modalities)
    
    vcf_files_1 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    )
    vcf_files_2 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells_2}/*_2.vcf.gz"
    )
    // Combine all VCF files
    vcf_files = vcf_files_1.mix(vcf_files_2)
        .collect()
        .view { files -> "All VCF files (${files.size()} total): ${files.collect { it.name }.join(', ')}" }

    index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_*/*.vcf.gz.csi"
    )
    
    // 1a. Filter variant calls
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        vcf_files,
        index_files.collect(),
        mod_channel.collect(),
        true
    )
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs, mod_channel.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs, mod_channel.collect())
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs, mod_channel.collect())
    HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())
}
#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variation on main_from_genotype.nf that can be used to assess
 * the impact of comparing pseudobulks to single-cell data on sample matching performance.
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
include { MERGE_AND_FILTER_VCFS as MERGE_PSEUDOBULK_VCFS} from './modules/1a_merge_vcfs'
include { MERGE_AND_FILTER_VCFS as MERGE_SC_VCFS} from './modules/1a_merge_vcfs'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'


// Workflow specific parameters
params.outdir = "results/pseudobulk_vs_sc"


// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)
    if (params.pseudobulk) {
        throw new IllegalArgumentException("main_pseudobulk_vs_sc.nf requires params.pseudobulk = false because it compares two data modes.")
    }
    requirePositiveInteger(params.ncells, 'params.ncells')

    validateVcfRunPaths(
        params,
        ['vireo', 'fingerprints', 'hysys', 'ngscheckmate'],
        [
            single_cell: "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*_single-cell.vcf.gz",
            pseudobulk: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
        ]
    )

    sc_vcf_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*_single-cell.vcf.gz"
    ).collect()
    .view { x -> "SC VCF Files: ${x}"}
    sc_index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*_single-cell.vcf.gz.csi"
    ).collect()

    pseudobulk_vcf_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    ).collect()
    .view { x -> "Pseudobulk VCF Files: ${x}"}
    pseudobulk_index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz.csi"
    ).collect()
    
    // 1a. Filter variant calls
    pseudobulk_channel = channel.value("pseudobulk")
    pseudobulk_modality_channel = pseudobulk_channel.map { x -> [x]}
    pseudobulk_channel.view { x -> "Pseudobulk channel: ${x}" }
    filtered_pseudobulk_vcfs = MERGE_PSEUDOBULK_VCFS(
        pseudobulk_vcf_files, 
        pseudobulk_index_files,
        pseudobulk_modality_channel,
        true
    )
    filtered_pseudobulk_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    sc_channel = channel.value("single-cell")
    sc_modality_channel = sc_channel.map { x -> [x]}
    sc_channel.view { x -> "Single-cell channel: ${x}" }
    filtered_sc_vcfs = MERGE_SC_VCFS(
        sc_vcf_files, 
        sc_index_files,
        sc_modality_channel,
        false
    )
    filtered_sc_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    modality_vcfs = filtered_pseudobulk_vcfs.modality_vcfs.mix(filtered_sc_vcfs.modality_vcfs).collect()
    individual_vcfs = filtered_pseudobulk_vcfs.individual_vcfs.mix(filtered_sc_vcfs.individual_vcfs).collect()
    modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    mod_channel = pseudobulk_modality_channel.mix(sc_modality_channel).collect()
    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(modality_vcfs, mod_channel)
    NGSCHECKMATE(individual_vcfs, mod_channel)
    CROSSCHECK_FINGERPRINTS(individual_vcfs, mod_channel)
    HYSYS(individual_vcfs, mod_channel)
}
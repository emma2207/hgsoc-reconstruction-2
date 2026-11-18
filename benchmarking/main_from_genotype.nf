#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variation on main.nf that can be used if variant calling 
 * has already been performed.
 * It analyzes VCF files to determine if samples come from the same donor
 * using multiple tools: HYSYS, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

include {
    isValidRealDataModalities;
    validateBenchmarkParameters
} from './modules/helpers/modality_utils'
include { validateVcfRunPaths } from './modules/helpers/path_validation'

// Import processes
include { MERGE_AND_FILTER_VCFS } from './modules/1a_merge_vcfs'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'
include { OMICSPRINT } from './modules/2a_omicsprint'
include { PEDDY } from './modules/2a_peddy'
include { TIMEATTACKGENCOMP } from './modules/2a_timeattackgencomp'


// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)

    if (!params.pseudobulk && !isValidRealDataModalities(params.modalities)) {
        throw new IllegalArgumentException(
            "Real-data runs require exactly two distinct, non-empty modality names in params.modalities, one for each modality to compare. Names must not contain one another."
        )
    }

    modalities = params.pseudobulk ? ["pseudobulk"] : params.modalities
    def vcfPattern = params.pseudobulk
        ? "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz"
        : "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*.vcf.gz"
    validateVcfRunPaths(
        params,
        ['vireo', 'fingerprints', 'hysys', 'ngscheckmate', 'peddy', 'timeattackgencomp', 'omicsprint'],
        ['sample VCF inputs': vcfPattern],
        params.pseudobulk ? [] : modalities
    )

    // Find VCFs for the selected experiment
    if (!params.pseudobulk) {
        vcf_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*.vcf.gz"
        )
        .view { x -> "VCF Files: ${x}"}
        index_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*.vcf.gz.csi"
        )
        .view { x -> "Index Files: ${x}"}
    } else {
        vcf_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz"
        )
        .view { x -> "VCF Files: ${x}"}
        index_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
        )
        .view { x -> "Index Files: ${x}"}
    }
    mod_channel = channel.fromList(modalities)
    
    // 1a. Filter variant calls
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        vcf_files.collect(), 
        index_files.collect(),
        mod_channel.collect(),
        params.pseudobulk
    )
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect(), modalities.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs, modalities.collect())
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs, modalities.collect())
    HYSYS(filtered_vcfs.individual_vcfs, modalities.collect())
    PEDDY(filtered_vcfs.all_variants, filtered_vcfs.all_variants_index, modalities.collect())
    TIMEATTACKGENCOMP(filtered_vcfs.individual_vcfs.collect(), filtered_vcfs.individual_vcfs_index.collect(), modalities)
    OMICSPRINT(filtered_vcfs.all_variants, modalities)
}
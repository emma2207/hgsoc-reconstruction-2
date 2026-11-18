#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variation on main_from_genotype.nf that can be used to assess
 * the impact of including samples multiple from the same donors multiple times 
 * on sample matching performance.
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
params.n_double_samples = 10  // Number of samples to randomly select and "duplicate" (different VCF files from the same donor)
params.repeat = 5  // Number of iterations to repeat the random sample duplication and analysis
params.outdir = "results/double_samples_${params.n_double_samples}/it_${params.repeat}"
  

// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)
    if (!params.pseudobulk) {
        throw new IllegalArgumentException("main_double_samples.nf requires params.pseudobulk = true.")
    }
    requirePositiveInteger(params.n_double_samples, 'params.n_double_samples')

    def inputFiles = validateVcfRunPaths(
        params,
        ['vireo', 'fingerprints', 'hysys', 'ngscheckmate'],
        [
            set_1: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz",
            set_2: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz",
            set_3: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_3.vcf.gz"
        ]
    )
    if (inputFiles.set_3.size() < params.n_double_samples) {
        throw new IllegalArgumentException(
            "params.n_double_samples (${params.n_double_samples}) exceeds the number of set-3 VCFs (${inputFiles.set_3.size()})."
        )
    }

    // Find all vcfs with the listed datatypes
    modalities = ["pseudobulk"]
    mod_channel = channel.fromList(modalities)
    
    // Use all pseudobulk files ending in _1.vcf.gz and _2.vcf.gz
    vcf_files_1 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    )
    vcf_files_2 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz"
    )
    vcf_files_3_all = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_3.vcf.gz"
    )
    .view { n -> "Number of VCF Files 3 before removing samples: ${n.size()}"}
    index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
    )

    // Select a random subset of vcf_files_3 to add as duplicate samples
    vcf_files_3 = vcf_files_3_all
    .collect()
    .map { files_3 ->
        def seededRandom = new Random(params.repeat + 100) // Seed inside closure for proper serialization
        def files_3_copy = files_3.collect()
        Collections.shuffle(files_3_copy, seededRandom)
        files_3_copy.take(params.n_double_samples)
    }
    .view { x -> "Number of VCF Files 3 added: ${x.size()}"}

    // Combine all VCF files
    vcf_files = vcf_files_1.mix(vcf_files_2).mix(vcf_files_3)
        .collect()
        .view { files -> "All VCF files (${files.size()} total): ${files.collect { it.name }.join(', ')}" }
    
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
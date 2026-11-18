#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variation on main_from_genotype.nf that can be used to assess
 * the impact of missing samples on sample matching performance.
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
params.n_samples_to_remove = 1  // Number of samples to randomly remove from one set of VCF files to simulate missing data
params.repeat = 5  // Number of iterations to repeat the random sample removal and analysis
params.outdir = "results/remove_samples_${params.n_samples_to_remove}/it_${params.repeat}"


// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)
    if (!params.pseudobulk) {
        throw new IllegalArgumentException("main_missing_samples.nf requires params.pseudobulk = true.")
    }
    requirePositiveInteger(params.n_samples_to_remove, 'params.n_samples_to_remove')

    def inputFiles = validateVcfRunPaths(
        params,
        ['vireo', 'fingerprints', 'hysys', 'ngscheckmate'],
        [
            set_1: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz",
            set_2: "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz"
        ]
    )
    if (inputFiles.set_2.size() <= params.n_samples_to_remove) {
        throw new IllegalArgumentException(
            "params.n_samples_to_remove (${params.n_samples_to_remove}) must be smaller than the number of set-2 VCFs (${inputFiles.set_2.size()}) so at least one sample remains."
        )
    }

    // Find all vcfs with the listed datatypes
    modalities = ["pseudobulk"]
    mod_channel = channel.fromList(modalities)
    
    // Use all pseudobulk files ending in _1.vcf.gz and _2.vcf.gz
    vcf_files_1 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    ).collect()
    vcf_files_2_all = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz"
    )
    .collect()
    .view { n -> "Number of VCF Files 2 before removing samples: ${n.size()}"}
    index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
    )

    // Remove random samples from vcf_files_2 to simulate missing data
    vcf_files_2 = vcf_files_2_all.map { files_2 ->
        def seededRandom = new Random(params.repeat) // Seed inside closure for proper serialization
        def files_2_copy = files_2.collect()
        Collections.shuffle(files_2_copy, seededRandom)
        files_2_copy.take(files_2_copy.size() - params.n_samples_to_remove)
    }
    .view { x -> "Number of VCF Files 2 after removing samples: ${x.size()}"}

    // Combine all VCF files
    vcf_files = vcf_files_1.mix(vcf_files_2)
        .flatMap { x -> x }
        .collect()
        .view { file -> "VCF file: ${file.name}" }
    
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
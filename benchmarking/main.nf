#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: BAMixChecker, HYSYS, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

include {
    isValidRealDataModalities;
    validateBenchmarkParameters
} from './modules/helpers/modality_utils'
include { validateBamRunPaths } from './modules/helpers/path_validation'

// Import processes
include { PREPARE_INPUTS } from './modules/0_prepare_inputs'
include { SNP_REFERENCE_INDEX } from './modules/0_snp_ref_index'
include { GENOTYPE_INDIVIDUAL } from './modules/1a_genotype_individual'
include { MERGE_AND_FILTER_VCFS } from './modules/1a_merge_vcfs'
include { FILTER_BAM } from './modules/1b_filter_bam'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'
include { PEDDY } from './modules/2a_peddy'
include { TIMEATTACKGENCOMP } from './modules/2a_timeattackgencomp'
include { OMICSPRINT } from './modules/2a_omicsprint'
include { BAMIXCHECKER } from './modules/2b_bamixchecker'
include { CONPAIR_PILEUP } from './modules/2b_conpair'
include { CONPAIR_VERIFY } from './modules/2b_conpair'
include { SOMALIER_EXTRACT } from './modules/2b_somalier'
include { SOMALIER_RELATE } from './modules/2b_somalier'


// Main workflow
workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)

    if (!params.pseudobulk && !isValidRealDataModalities(params.modalities)) {
        throw new IllegalArgumentException(
            "Real-data runs require exactly two distinct, non-empty modality names in params.modalities, one for each modality to compare. Names must not contain one another."
        )
    }

    modalities = params.pseudobulk ? ["pseudobulk"] : params.modalities
    validateBamRunPaths(params, [
        'vireo',
        'fingerprints',
        'hysys',
        'ngscheckmate',
        'peddy',
        'timeattackgencomp',
        'omicsprint',
        'bamixchecker',
        'conpair',
        'somalier'
    ])

    // Find all bams with the listed datatypes
    if (!params.pseudobulk) {
        bam_patterns = modalities.collect { datatype ->
            "${params.bamsDir}/${params.dataset}/${datatype}/*/Aligned.sortedByCoord.out.bam" 
        }
        
        bam_files = channel
            .fromPath(bam_patterns)
            .map { f ->
                // safe: get the last element of the parent Path as String
                def parentName = f.parent ? f.parent.getFileName().toString() : ''
                def datatype = modalities.find { modality -> f.toString().contains("/${modality}/") } ?: 'unknown'
                tuple(parentName, f.baseName, f, datatype)
            }
            .view { x -> "Found BAM files: ${x}" }
    } else {
        bam_files = channel.fromPath("${params.pseudobulkDir}/${params.dataset}/*/pseudobulk_n_barcodes_${params.ncells}_*.bam")
            .map { f ->
                // safe: get the last element of the parent Path as String
                def parentName = f.parent ? f.parent.getFileName().toString() : ''
                tuple(parentName, f.baseName, f, 'pseudobulk')
            }
            .view { x -> "Found Pseudobulk BAM files: ${x}" }
    }
    
    // 0. Prepare inputs
    prepped_bam_files = PREPARE_INPUTS(bam_files)
    prepped_bam_files.view { x -> "Prepared BAM files: ${x}" }
    // 0. Prepare SNP reference
    snp_ref = SNP_REFERENCE_INDEX(params.snp_vcf)

    // 1a. Genotype individual samples in parallel
    individual_vcfs = GENOTYPE_INDIVIDUAL(prepped_bam_files, snp_ref.vcf, snp_ref.tbi)

    // 1a. Merge and filter variants
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        individual_vcfs.vcf.collect(), 
        individual_vcfs.index.collect(), 
        modalities,
        params.pseudobulk
    )
    filtered_vcfs.modality_vcfs.view { x -> "Modality VCFs: ${x}"}
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 1b. Filter aligned reads
    filtered_bams = FILTER_BAM(bam_files)
    filtered_bams.bam.collect().view { x -> "Filtered BAMs: ${x}" }

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect(), modalities.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs, modalities.collect())
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs, modalities.collect())
    HYSYS(filtered_vcfs.individual_vcfs, modalities.collect())
    PEDDY(filtered_vcfs.all_variants, filtered_vcfs.all_variants_index, modalities.collect())
    TIMEATTACKGENCOMP(filtered_vcfs.individual_vcfs.collect(), filtered_vcfs.individual_vcfs_index.collect(), modalities)
    OMICSPRINT(filtered_vcfs.all_variants, modalities)

    // 2b. Run similarity analysis tools in parallel on filtered BAMs
    BAMIXCHECKER(filtered_bams.bam.collect(), modalities)
    conpair_pairs = filtered_bams.bam.collect(flat: false).flatMap { bams ->
            if (params.pseudobulk) {
                def bams_1 = bams.findAll { bam -> bam.baseName.toString() ==~ /.*_1_filtered$/ }
                def bams_2 = bams.findAll { bam -> bam.baseName.toString() ==~ /.*_2_filtered$/ }

                if (bams_1.isEmpty() || bams_2.isEmpty()) {
                    throw new IllegalStateException("CONPAIR pseudobulk pairing failed: expected both *_1_filtered.bam and *_2_filtered.bam inputs")
                }

                bams_1.withIndex().collectMany { bam_1, i ->
                    bams_2.withIndex().collect { bam_2, j ->
                        def pair_id = "pb_${i}_${j}"
                        def sample_1 = "${pair_id}_A_${bam_1.baseName}"
                        def sample_2 = "${pair_id}_B_${bam_2.baseName}"
                        tuple(sample_1.toString(), bam_1, sample_2.toString(), bam_2, "pseudobulk_1", "pseudobulk_2")
                    }
                }
            } else {
                def modality_1 = modalities[0]
                def modality_2 = modalities[1]
                def bams_mod_1 = bams.findAll { bam -> bam.baseName.toString().contains("_${modality_1}_filtered") }
                def bams_mod_2 = bams.findAll { bam -> bam.baseName.toString().contains("_${modality_2}_filtered") }

                if (bams_mod_1.isEmpty() || bams_mod_2.isEmpty()) {
                    throw new IllegalStateException("CONPAIR real-data pairing failed: expected BAMs from both modalities '${modality_1}' and '${modality_2}'")
                }

                bams_mod_1.withIndex().collectMany { bam_1, i ->
                    bams_mod_2.withIndex().collect { bam_2, j ->
                        def pair_id = "rd_${i}_${j}"
                        def sample_1 = "${pair_id}_A_${bam_1.baseName}"
                        def sample_2 = "${pair_id}_B_${bam_2.baseName}"
                        tuple(sample_1.toString(), bam_1, sample_2.toString(), bam_2, modality_1, modality_2)
                    }
                }
            }
    }
    conpair_all_pairs = conpair_pairs.collect(flat: false)
    conpair_pileup_inputs = conpair_all_pairs.flatMap { all_pairs ->
        def bam_by_path = [:]
        all_pairs.each { pair ->
            bam_by_path[pair[1].toString()] = pair[1]
            bam_by_path[pair[3].toString()] = pair[3]
        }

        bam_by_path.entrySet().toList().withIndex().collect { entry, index ->
            tuple(entry.key, "bam_${index}", entry.value)
        }
    }

    CONPAIR_PILEUP(conpair_pileup_inputs)

    conpair_pileup_map = CONPAIR_PILEUP.out.pileup
        .collect(flat: false)
        .map { records -> records.collectEntries { record -> [(record[0]): record[1]] } }

    conpair_pair_contexts = conpair_all_pairs.flatMap { all_pairs ->
        all_pairs.collect { pair ->
            tuple(
                pair[1].toString(), pair[0],
                pair[3].toString(), pair[2],
                pair[4], pair[5]
            )
        }
    }

    conpair_verification_inputs = conpair_pair_contexts
        .combine(conpair_pileup_map)
        .map { bam_key_1, sample_1, bam_key_2, sample_2, mod_1, mod_2, pileups_by_key ->
            tuple(
                sample_1, pileups_by_key[bam_key_1],
                sample_2, pileups_by_key[bam_key_2],
                mod_1, mod_2
            )
        }

    CONPAIR_VERIFY(conpair_verification_inputs)
    somalier_files = SOMALIER_EXTRACT(filtered_bams.bam, filtered_bams.bam_index, modalities)
    SOMALIER_RELATE(somalier_files.collect(), modalities[0], modalities.size() > 1 ? modalities[1] : modalities[0])
}
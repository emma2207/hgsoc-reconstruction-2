# Benchmarking Nextflow Pipeline

Nextflow pipeline to test tools to detect sample mislabeling (i.e. sample swaps or match samples).
This workflow is branched into two paths: one for .vcf files getting fed into tools and one for tools that need .bam files.

## Main files
- **sample-matching.yml:** yml specifying software dependencies for everything but the individual tools.
- **run_nf_pipeline.sh:** shell script to run the whole pipeline on the cluster
- **nextflow.config:** specifies global parameters and profiles for different execution environments. Note that the params can easily be overwritten via a command line argument at the time of execution.
- **main.nf:** main nextflow file. This is the file that needs to be executed to run the whole workflow. This particular main.nf is the most basic one, starting from .bam files.
    - **main_from_genotype.nf:** a variation on `main.nf` where the input is .vcf instead of .bam. I.e. the genotyping has already been performed.
    - **main_missing_samples.nf:** a variation on `main_from_genotype.nf` for pseudobulk experiments where the input data is manipulated to have missing samples.
    - **main_double_samples.nf:** a variation on `main_from_genotype.nf` for pseudobulk experiments where the input data is manipulated to have samples that are doubly represented (not literally identical samples, but drawn from the same underlying sample).
    - **main_uneven_pseudobulk_sizes.nf:** variation on `main_from_genotype.nf` for pseudobulk experiments where pseudobulk drawn from the size samples but with different numbers of cells are compared.
    - **main_pseudobulk_vs_sc.nf:** variation on `main_from_genotype.nf` comparing pseudobulks with single-cell (or single-nucleus) samples.

## Software versions
The Conda YAML files pin the versions of their declared dependencies to those recorded in the HPC environments used for the analysis. 

The pipeline currently points processes to existing environments under `params.conda`; it does not automatically create them from these YAML files.

## Modules
Each module specifies a nextflow process that is called in main:
- **0_prepare_inputs.nf:** prepare input files .bam files
- **0_snp_ref_index.nf:** prepares the SNP reference index.
- **1a_genotype_individual.nf:** performs variant calling for each sample. It results in a .vcf file and these get fed into modules whose file starts with 2a.
- **1a_filter_vcf.nf:** filters VCF inputs.
- **1a_merge_vcfs.nf:** merges and filters individual and modality VCFs.
- **1b_filter_bam.nf:** filters the aligned reads by quality score and these filtered .bam's get fed into any module whose file starts with 2b.
- **2a_fingerprints.nf:** applies the tool CrosscheckFingerprints and needs the conda env specified in `fingerprints.yml`.
- **2a_hysys.nf:** applies the tool HaveYouSwappedYourSamples and needs the conda env specified in `hysys.yml`.
- **2a_ngscheckmate.nf:** applies the tool NGSCheckmate and needs the conda env specified in `ngscheckmate.yml`.
- **2a_vireo.nf:** applies the tool Vireo through a python script `vireo.py`.
- **2a_omicsprint.nf:** applies the Bioconductor package omicsPrint in R on SNPs extracted from the merged VCF.
- **2a_peddy.nf:** applies peddy to VCF data using `peddy.yml`.
- **2a_timeattackgencomp.nf:** applies TimeAttackGenComp to VCF data using `timeattackgencomp.yml`.
- **2b_bamixchecker.nf:** applies BAMixChecker and uses the conda env specified in `bamixchecker.yml`.
- **2b_conpair.nf:** runs `CONPAIR_PILEUP` once per unique BAM, then `CONPAIR_VERIFY` for each pairwise comparison, using `conpair.yml`.
- **2b_somalier.nf:** runs Somalier extraction and pairwise relatedness analysis using `somalier.yml`.

## Parameters
Some notes on the parameter that need to be set in `nextflow.config`.
Input and reference paths are set with `params.bamsDir`, `params.pseudobulkDir`, `params.vcfsDir`, `params.refGenome`, and `params.snp_vcf`; `params.outdir` sets the output location. Adjust these paths for the environment as needed.

### Pseudobulk experiments
- **params.dataset:** name of dataset; defaults to `hgsoc-new`. The config has dataset-specific read-depth entries for `hgsoc`, `hgsoc-new`, `high_grade_glioma`, `low_grade_glioma`, and `wilms_tumor`; other labels use the fallback settings.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `true` for pseudobulk experiments; the config default is `false`.
- **params.ncells:** number of cells included in the pseudobulks; defaults to `50000` and must match the input filenames. Common dataset sizes include 10,000; 50,000; 100,000; and 500,000 cells.
- **params.read_depth:** read-depth thresholds are configured per dataset and modality, with an `upper` cap and fallback. See the values in `nextflow.config`.
- **params.mpileup_max_depth:** maximum mpileup depth for genotype calling; defaults to `100`.
- **params.omicsprint_max_snps:** optional cap on number of SNPs passed to omicsPrint after filtering. Use `0` to keep all SNPs.
- **params.omicsprint_call_rate:** SNP call-rate threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.80`; omicsPrint default: `0.95`). Lower values keep more SNPs in sparse datasets.
- **params.omicsprint_coverage_rate:** sample coverage threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.25`; omicsPrint default: `2/3`). Lower values retain more samples in low-coverage comparisons.
- **params.timeattackgencomp_target_bed:** target BED used by TimeAttackGenComp.
- **Conpair parameters:** `params.conpair_marker_bed`, `params.conpair_marker_txt`, `params.conpair_min_cov`, `params.conpair_min_map_qual`, and `params.conpair_min_base_qual` override marker files and concordance filters.

### Real data experiments
- **params.dataset:** name of dataset; defaults to `hgsoc-new`. The config has dataset-specific read-depth entries for `hgsoc`, `hgsoc-new`, `high_grade_glioma`, `low_grade_glioma`, and `wilms_tumor`; other labels use the fallback settings.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `false` for real data experiments; this is the config default.
- **params.ncells:** number of cells included in the pseudobulks. This parameter is not used when `params.pseudobulk = false`.
- **params.read_depth:** read-depth thresholds are configured per dataset and modality, with an `upper` cap and fallback. See the values in `nextflow.config`.
- **params.mpileup_max_depth:** maximum mpileup depth for genotype calling; defaults to `100`.
- **params.omicsprint_max_snps:** optional cap on number of SNPs passed to omicsPrint after filtering. Use `0` to keep all SNPs.
- **params.omicsprint_call_rate:** SNP call-rate threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.80`; omicsPrint default: `0.95`). Lower values keep more SNPs in sparse datasets.
- **params.omicsprint_coverage_rate:** sample coverage threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.25`; omicsPrint default: `2/3`). Lower values retain more samples in low-coverage comparisons.
- **params.timeattackgencomp_target_bed:** target BED used by TimeAttackGenComp.
- **Conpair parameters:** `params.conpair_marker_bed`, `params.conpair_marker_txt`, `params.conpair_min_cov`, `params.conpair_min_map_qual`, and `params.conpair_min_base_qual` override marker files and concordance filters.

Additionally, for real data experiments, one has to provide the two modalities that one wants to compare. 
This is set in whichever main Nextflow files one is using (look for `modalities`).
The options are specific to each dataset:
- **hgsoc:** bulk_chunk_ribo, bulk_dissociated_polyA, bulk_dissociated_ribo & single-cell.
- **hgsoc-new:** bulk_chunk_ribo & bulk_diss_polyA.
- **high_grade_glioma:** bulk & single-cell.
- **low_grade_glioma:** bulk & single-cell.
- **wilms_tumor:** bulk & single-nucleus.

## ntsm
The `ntsm/` folder contains a standalone Nextflow workflow for running NTSM. It maps real-data samples to single- or paired-end FASTQs from SRA metadata, or converts pseudobulk BAMs to FASTQ, then runs NTSM counts and pairwise evaluation. Its entry point, settings, and Conda environment are in `ntsm/main.nf`, `ntsm/nextflow.config`, and `ntsm/ntsm.yml`; launch it with `ntsm/run_nf_pipeline.sh`.

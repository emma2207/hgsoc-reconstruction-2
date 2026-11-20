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
The original pipeline run used Nextflow **24.04.4**. The Conda YAML files
pin the versions of their declared dependencies to those recorded in the HPC
environments used for the analysis.

The pipeline uses existing environments under `params.conda`; it does not
automatically create them from these YAML files. Create the environments
yourself and set `params.conda` and the module `conda` directives to match
their locations.

## Using your own dataset
See [Using your own dataset](USING_YOUR_OWN_DATA.md) for workflow selection,
input layouts, reference resources, sample naming, and dataset-specific
settings.

## Running the pipeline
Run the commands below from the `benchmarking/` directory so that Nextflow loads this directory's `nextflow.config`. 
First install or load Nextflow and Conda, create the required Conda environments, and adjust the input reference, tool, environment, and output paths in `nextflow.config` for your system. 
See **Parameters** below for dataset settings. 
The existing paths and scheduler settings are specific to the original HPC environment.
Each workflow runs preflight checks for its required inputs, indexes, tool
resources, and writable output directory before starting processes. Conda
environment locations and contents are not checked; see
[Using your own dataset](USING_YOUR_OWN_DATA.md#preflight-validation).

### Using the launch script on Slurm
Before submitting `run_nf_pipeline.sh`:

1. Change the `#SBATCH` account, partition, QoS, time, memory, and email settings to match your cluster, or remove directives your cluster does not use. 
These settings apply to the job running Nextflow itself.
2. In the `slurm` profile in `nextflow.config`, adjust `process.queue`, `process.clusterOptions`, and resource settings for the jobs Nextflow submits. 
These are separate from the script's `#SBATCH` settings.
3. Adapt the `module load` commands to the modules on your cluster. Review the Alpine-specific `PATH`, `CPATH`, and `C_INCLUDE_PATH` exports and remove or change them if they are not appropriate for your system.

Then submit the script with the desired workflow file as an argument:

```bash
sbatch run_nf_pipeline.sh main.nf
```

The script uses the submission directory as its working directory on Slurm, so submit it from `benchmarking/`. 
It selects the `slurm` profile, uses `benchmarking/nextflow/` for the work directory, enables task output and a trace report, and passes `-resume` to reuse eligible cached tasks.

### Running Nextflow directly
Alternatively, after loading Nextflow and Conda and adapting the Slurm profile as above, launch the workflow directly:

```bash
nextflow run main.nf -profile slurm
```

This runs the Nextflow controller in your current session and submits process jobs to Slurm. 
Follow your cluster's policy for long-running controllers; you may need an allocation or a batch script rather than a login node. 
To match the launch script's work directory, logging, and resume options, use:

```bash
nextflow run main.nf -profile slurm -w nextflow -process.echo -resume -with-trace
```

Replace `main.nf` with another workflow file listed above when appropriate.
For local execution, use `-profile local` instead, after adjusting the local profile's paths and resource settings. 
Running `bash run_nf_pipeline.sh main.nf` locally also selects that profile, but its unconditional `module load` commands must first be adapted or removed if your machine does not provide those modules.

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

## Expected output files
The modules below publish their declared outputs under the corresponding
subdirectory of `params.outdir`, usually retaining a nested
`${params.dataset}/...` directory structure. Names in square brackets are
generated from sample IDs, modalities, or read-depth settings. Outputs marked
optional may not be produced for every run. Intermediate files not declared
as process outputs (including task logs and temporary files) remain in the
Nextflow work directory.

| Module / process | Published files |
| --- | --- |
| `0_prepare_inputs.nf` / `PREPARE_INPUTS` | Renamed BAMs and their `.bam.bai` indexes under `0_prepared_bams/`. Names follow `[sample_id]_[nr].bam`, where `nr` is supplied by the workflow input tuple. |
| `0_snp_ref_index.nf` / `SNP_REFERENCE_INDEX` | `common_snps.vcf.gz`, `common_snps.vcf.gz.tbi`, and `versions.yml` under `0_snp_ref_index/`. |
| `1a_genotype_individual.nf` / `GENOTYPE_INDIVIDUAL` | `[sample].vcf.gz` and `[sample].vcf.gz.csi` under `1a_individual_vcf/`. |
| `1a_filter_vcf.nf` / `FILTER_VCF` | Per-sample `*_individual_variants.vcf.gz`; per-modality `*_modality_variants.vcf.gz`; and `filtered_variants_[modality]_rd_[threshold].vcf.gz` under `1a_vcf/`. |
| `1a_merge_vcfs.nf` / `MERGE_AND_FILTER_VCFS` | `all_variants.vcf.gz` and its `.csi` index; per-sample `*_individual_variants.vcf.gz` and `.csi` indexes; per-modality `*_modality_variants.vcf.gz`; and `filtered_variants_[modality]_rd_[threshold].vcf.gz` under `1a_vcf/`. |
| `1b_filter_bam.nf` / `FILTER_BAM` | `[sample]_filtered.bam` and `[sample]_filtered.bam.bai` under `1b_filtered_bam/`. |
| `2a_fingerprints.nf` / `CROSSCHECK_FINGERPRINTS` | `crosscheck_metrics.txt` under `2a_fingerprints/`. |
| `2a_hysys.nf` / `HYSYS` | `concordance_output.txt` and `model_results.txt` under `2a_hysys/`. |
| `2a_ngscheckmate.nf` / `NGSCHECKMATE` | `output_all.txt`, `output_matched.txt`, and `output_output_corr_matrix.txt` under `2a_ngscheckmate/`. |
| `2a_vireo.nf` / `VIREO_MATCH` | `similarity_matrix.csv` and `matched_samples.csv` under `2a_vireo/`. |
| `2a_omicsprint.nf` / `OMICSPRINT` | `omicsprint_genotype_matrix.tsv`, `omicsprint_input_samples.txt`, `omicsprint_allele_sharing.tsv`, `omicsprint_mismatches.tsv`, and `omicsprint_ibs_scatter.pdf` under `2a_omicsprint/`. |
| `2a_peddy.nf` / `PEDDY` | Files matching `peddy*` under `2a_peddy/`. |
| `2a_timeattackgencomp.nf` / `TIMEATTACKGENCOMP` | `timeattackgencomp.snv.out.txt`, `timeattackgencomp.raw.out`, and optionally `timeattackgencomp.pdf` under `2a_timeattackgencomp/`. |
| `2b_bamixchecker.nf` / `BAMIXCHECKER` | `BAMixChecker_Report.html`, `BAMixChecker_Heatmap.pdf`, and `Total_result.txt`; `Matched_samples.txt` and `Mismatched_samples.txt` are optional, under `2b_bamixchecker/`. |
| `2b_conpair.nf` / `CONPAIR_PILEUP` | `[sample]_pileup.txt` is an intermediate passed to `CONPAIR_VERIFY`; it is not published to `params.outdir`. |
| `2b_conpair.nf` / `CONPAIR_VERIFY` | `[sample_1]_vs_[sample_2]_concordance.txt` under `2b_conpair/`. |
| `2b_somalier.nf` / `SOMALIER_EXTRACT` | `[sample].somalier` files under `2b_somalier/extract/`. |
| `2b_somalier.nf` / `SOMALIER_RELATE` | `somalier*.tsv` and optionally `somalier*.html` under `2b_somalier/`. |

Not every module runs in every workflow: for example, BAM-based tools are not
used by the existing-VCF entry point. The separate NTSM workflow has its own
outputs and configuration under `ntsm/`.

## Parameters
Some notes on the parameter that need to be set in `nextflow.config`.
Input and reference paths are set with `params.bamsDir`, `params.pseudobulkDir`, `params.vcfsDir`, `params.refGenome`, and `params.snp_vcf`; `params.outdir` sets the output location. Adjust these paths for the environment as needed.

### Shared experiment parameters
- **`params.dataset`:** dataset label; defaults to `hgsoc-new`. Dataset-specific read-depth entries are defined for `hgsoc`, `hgsoc-new`, `high_grade_glioma`, `low_grade_glioma`, and `wilms_tumor`. Other labels use the `default` read-depth entry.
- **`params.pseudobulk`:** selects pseudobulk (`true`) or real-data (`false`) inputs; defaults to `false`.
- **`params.read_depth`:** read-depth thresholds are configured per dataset and modality, with an `upper` cap and fallback. Review and set values appropriate for your data in `nextflow.config`.
- **`params.mpileup_max_depth`:** maximum mpileup depth for genotype calling; defaults to `100`.
- **`params.omicsprint_max_snps`:** optional cap on SNPs passed to omicsPrint after filtering. Use `0` to keep all SNPs.
- **`params.omicsprint_call_rate`:** SNP call-rate threshold passed to `omicsPrint::alleleSharing` (pipeline default `0.80`; omicsPrint default `0.95`). Lower values retain more SNPs in sparse datasets.
- **`params.omicsprint_coverage_rate`:** sample coverage threshold passed to `omicsPrint::alleleSharing` (pipeline default `0.25`; omicsPrint default `2/3`). Lower values retain more samples in low-coverage comparisons.
- **`params.timeattackgencomp_target_bed`:** target BED used by TimeAttackGenComp.
- **Conpair settings:** `params.conpair_marker_bed`, `params.conpair_marker_txt`, `params.conpair_min_cov`, `params.conpair_min_map_qual`, and `params.conpair_min_base_qual` override marker files and concordance filters.

### Mode-specific parameters
- **Pseudobulk:** set **`params.ncells`** to the number of cells in the pseudobulks (default `50000`); it must match the input directory and filenames. Common sizes include 10,000; 50,000; 100,000; and 500,000 cells. In the pseudobulk-specific experiment workflows, the modality label is fixed to `pseudobulk`; `params.modalities` does not replace it.
- **Real data:** `params.ncells` is not used. For `main.nf` and `main_from_genotype.nf`, set **`params.modalities`** to exactly two distinct, non-empty modality names in `nextflow.config`. Both names must match the input layout and VCF/BAM sample naming. Other experiment workflows have fixed modalities; see [Using your own dataset](USING_YOUR_OWN_DATA.md).

The configured modality names and dataset read-depth entries are not fully consistent: for `hgsoc-new`, the configured name is `bulk_diss_polyA`, while the read-depth map uses `bulk_dissociated_polyA`. If you change the read-depth threshold for this modality, make the map key match the name in `params.modalities`. Existing values for both are `0`, so the mismatch currently leaves the effective threshold at zero.

Modality labels used by the existing datasets include:
- **`hgsoc`:** `bulk_chunk_ribo`, `bulk_dissociated_polyA`, `bulk_dissociated_ribo`, `single-cell`.
- **`hgsoc-new`:** `bulk_chunk_ribo`, `bulk_diss_polyA`.
- **`high_grade_glioma` / `low_grade_glioma`:** `bulk`, `single-cell`.
- **`wilms_tumor`:** `bulk`, `single-nucleus`.

## ntsm
The `ntsm/` folder contains a standalone Nextflow workflow for running NTSM. 
It maps real-data samples to single- or paired-end FASTQs from SRA metadata, or converts pseudobulk BAMs to FASTQ, then runs NTSM counts and pairwise evaluation. 
Its entry point, settings, and Conda environment are in `ntsm/main.nf`, `ntsm/nextflow.config`, and `ntsm/ntsm.yml`; launch it with `ntsm/run_nf_pipeline.sh`.
The original NTSM pipeline run also used Nextflow **24.04.4**. Its
`hgsoc-new` branch selects the organized FASTQ input layout; it does not
exclude samples by ID.
For real-data runs, set `params.modalities` in `ntsm/nextflow.config` to
the two modality directory/metadata labels for the dataset. The
metadata-driven branch uses these names to locate its metadata CSVs and FASTQ
directories; the `hgsoc-new` branch uses them as labels for NTSM evaluation.
For `hgsoc-new`, set them to the names of the modality directories under
`params.pl_hgsoc`, for example `["bulk_chunk_ribo", "bulk_diss_polyA"]`.
Pseudobulk runs continue to use the fixed `pseudobulk` label.

### NTSM preflight validation
Before launching processes, NTSM checks that `params.dataset` is a safe
directory name, `params.pseudobulk` is a boolean, `params.threads` is a
positive integer, and `params.read_type` is `single` or `paired`. Pseudobulk
requires a positive integer `params.ncells`; real data requires two distinct,
safe modality names, neither containing the other, and a safe `params.datatype`
directory name.

The workflow checks the sites FASTA and creates and tests `params.outdir`.
For pseudobulk, it checks that the configured BAM pattern matches non-empty,
readable files. For the `hgsoc-new` layout, it checks the R1 and R2 files,
their pairing, and their modality directories; this branch requires both
mates even when `params.read_type` is `single`. For other real data, it
checks each metadata CSV, FASTQ directory, and metadata helper script.
After metadata processing, it checks that each modality has sample mappings,
that sample names use only letters,
numbers, `.`, `_`, or `-`, and that each mapped run has readable, non-empty
FASTQs, including R2 for paired-end data.

These checks do not validate file contents, sample identity, or reference
compatibility. As with the main workflow, Conda environment paths and
executables are not checked; configure the environment yourself.

### Expected NTSM output files
Published outputs are placed under `params.outdir` (default: `results`):

| Module / process | Published files |
| --- | --- |
| `modules/process_metadata.nf` / `PROCESS_METADATA` | `metadata/[dataset]/[read_type]/run_name_mapping_[modality].csv`, associating samples with FASTQ runs. Separate dataset, read-type, and modality names prevent mappings from overwriting one another. |
| `modules/pseudobulk_bam_to_fastq.nf` / `PSEUDOBULK_BAM_TO_FASTQ` | `[sample].R1.fastq.gz` and `[sample].R2.fastq.gz` under `pseudobulk_fastq/`. In single-end mode, R2 is an empty placeholder. |
| `modules/ntsm_count.nf` / `NTSM_COUNT` | `counts_[modality]_[sample].txt` under `ntsm_counts/`. |
| `modules/ntsm_eval.nf` / `NTSM_EVAL` | `ntsm_pairwise_pseudobulk_[dataset]_ncells_[ncells].tsv` for pseudobulk, or `ntsm_pairwise_[modalities]_[dataset].tsv` for real data, under `ntsm_eval/`. |

`NTSM_COUNT` creates concatenated per-sample FASTQs in its task directory,
but they are not declared outputs and are not published as final results.
The metadata mapping is only generated for the metadata-driven real-data
branch; pseudobulk and `hgsoc-new` use their configured BAM or organized
FASTQ inputs instead.

For direct execution, run from `benchmarking/ntsm/` so that its own `nextflow.config` is loaded:

```bash
nextflow run main.nf -profile slurm
```

Adapt that config's paths and Slurm settings independently of the main benchmarking config. 
The NTSM launch script has its own `#SBATCH` directives and module settings to adapt.
When using `sbatch`, submit from `benchmarking/ntsm/`; the script uses Slurm's
`SLURM_SUBMIT_DIR` as its working directory.
For local execution, its default executor is local and no profile is needed.

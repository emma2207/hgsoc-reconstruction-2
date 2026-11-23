# Using your own dataset

This guide covers the main benchmarking workflows in `benchmarking/`. They
were written for the specific datasets used in this repository and are not yet generic
sample-sheet-driven pipelines: using another dataset can require changing
workflow code as well as configuration. The examples below illustrate the
expected inputs; they are not a tested example dataset. The standalone NTSM
workflow has a separate configuration and is not covered here.

Run commands from `benchmarking/`. Before starting, follow the
[running instructions](README.md#running-the-pipeline) to set up Nextflow,
Conda, software and the execution profile. Check the pipeline's output using
non-sensitive test data before running a full dataset.
The original pipeline run used Nextflow **24.04.4**; see the README's
[software version notes](README.md#software-versions).

## 1. Choose a workflow

| Input and experiment | Entry point | Tools |
| --- | --- | --- |
| Coordinate-sorted BAMs, real data or pseudobulk | `main.nf` | VCF- and BAM-based tools |
| Individual VCFs, real data or pseudobulk | `main_from_genotype.nf` | VCF-based tools only |

Start with `main.nf` if genotyping is needed, or
`main_from_genotype.nf` if you already have individual VCFs. The other
entry points implement specific experiments and have additional filename
and parameter assumptions described below.

## 2. Configure paths and reference resources

In [`nextflow.config`](nextflow.config), set:

- `params.dataset`: a short dataset label used as a directory name.
- `params.pseudobulk`: `false` for real-data modalities or `true` for
  pseudobulk. For pseudobulk also set `params.ncells` to the cell count in
  the filenames; it must be a positive integer. The pipeline requires a
  non-empty dataset label and a boolean `params.pseudobulk`.
- For real-data runs of `main.nf` or `main_from_genotype.nf`,
  `params.modalities` must contain exactly two distinct, non-empty names.
  Since tools group samples by substring matching, neither name may contain
  the other.
- `params.bamsDir` or `params.pseudobulkDir` for BAM input; `params.vcfsDir`
  for existing VCF input. These are **parent directories**: the workflow
  appends the dataset and other subdirectories shown below.
- `params.outdir`: a writable output directory. Leave `params.projectDir`
  pointing to the benchmarking directory, where modules and helper scripts
  are located.
- `params.conda`: set this to the parent directory containing the Conda
  environments. The workflow expects specific environment subdirectory
  names, so keep those names or update the `conda` directives in the
  relevant modules. The preflight does not validate Conda environment
  locations or contents. The tools below were obtained from GitHub repositories;
  set their paths to the corresponding cloned directories, not just their
  Conda dependencies:

  | Parameter | GitHub repository |
  | --- | --- |
  | `params.CONPAIR` | [Conpair](https://github.com/nygenome/Conpair) |
  | `params.HYSYS` | [HaveYouSwappedYourSamples](https://github.com/PapenfussLab/HaveYouSwappedYourSamples) |
  | `params.NGS_CHECKMATE` | [NGSCheckMate](https://github.com/parklab/NGSCheckMate) |
  | `params.BAMIXCHECKER` | [BAMixChecker](https://github.com/heinc1010/BAMixChecker) |
  | `params.TIMEATTACKGENCOMP` | [TimeAttackGenComp](https://github.com/teerjk/TimeAttackGenComp) |

  For example, clone Conpair with
  `git clone https://github.com/nygenome/Conpair.git /path/to/Conpair`
  and set `params.CONPAIR` to `/path/to/Conpair`. Apply the same approach
  to the other repositories, following their installation instructions.
  `params.NGS_CHECKMATE` must point to the directory containing `ncm.py`
  and `SNP/`; its default is inside the `ngscheckmate` environment, but a
  clone elsewhere can be used by changing the parameter.

  Record the commit used for each clone, as the YAML pins do not fix the
  repository versions. **We do not recommend attempting to run BAMixChecker:**
  the version used for this analysis required a manual upgrade to Python 3
  to run. The upstream repository documents Python 2.7, and a fresh clone
  will not necessarily work with the pinned environment. `main.nf` currently
  calls BAMixChecker, so running the other tools without it requires
  disabling that process call in the workflow.

For BAM input, the selected profile sets `params.refGenome` to a **directory**
containing `fasta/genome.fa`, not to the FASTA itself. Provide an indexed
FASTA and the reference files required by the chosen tools, such as a
sequence dictionary. Set `params.snp_vcf` to a common-SNP VCF compatible
with that reference. Both settings are defined in the `local` and `slurm`
profiles, so update the profile you use or override them on the command
line.

For the original analysis, `refdata-gex-GRCh38-2024-A` was downloaded from
the [Cell Ranger reference website](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes).

Use the same genome assembly and compatible chromosome names in BAMs, VCFs,
FASTA, SNP and marker files. The SNP-reference preparation step renames
numbered chromosomes to `chr` names; check that this matches your reference.
The existing tool resources are GRCh38-oriented: in particular, review the
NGSCheckMate SNP BED, Conpair marker BED/TXT,
`params.timeattackgencomp_target_bed`, and the Picard haplotype map at
`modules/hg38_chr.map`. The Picard module expects that exact path; supply a
compatible map there or change the module. These resources may require
different files or code changes for another assembly.

## Preflight validation

Each workflow checks its required paths and inputs before starting processes.
It reports unmatched input patterns, unreadable or empty input files,
missing CSI indexes for individual input VCFs, and missing required
reference/tool files. The SNP reference VCF in `params.snp_vcf` does not
need a pre-existing index; the pipeline prepares an indexed copy.
It does not check Conda environment locations, contents, or executables;
create the environments yourself and ensure the module `conda` directives
point to them.
For BAM workflows, the pipeline indexes input BAMs itself, so a pre-existing
BAM index is not required. The preflight also creates `params.outdir` if
needed and tests that it can write there.

`main.nf` currently invokes all of its tools, including BAMixChecker, so its
preflight requires all of their environments and resources. If you disable a
tool call in `main.nf`, also remove that tool from the `validateBamRunPaths`
tool list or its environment and resources will still be required.

These checks do not confirm that BAM/VCF contents are valid, that sample IDs
are unique or correctly labeled, or that reference builds and contig names
are compatible. Check those properties separately; see **Check assumptions
and run** below.

## 3. Arrange the inputs

There is currently no sample sheet for these workflows. File paths and
names determine which samples are discovered and how they are grouped.
Use distinct, shell-safe sample names without spaces. Avoid modality names
that are substrings of each other or of unrelated sample names: several
processes group samples by substring matching.

### BAM input: `main.nf`

For real data, set `params.modalities` in `nextflow.config` to the **two**
modality names you want to compare. Both `main.nf` and
`main_from_genotype.nf` use this pair for real data. For example, with
`params.dataset = "my_dataset"` and
`params.modalities = ["bulk", "single-cell"]`, the directory layout is:

```text
<bamsDir>/
  my_dataset/
    bulk/
      sampleA/Aligned.sortedByCoord.out.bam
    single-cell/
      sampleA/Aligned.sortedByCoord.out.bam
```

There can be multiple sample directories per modality. Inputs must be
coordinate-sorted BAMs; the pipeline indexes them. It derives the sample
name from the parent directory and the modality from the path. If your BAMs
use different names or directory layouts, update both the discovery globs in
`main.nf` and the matching patterns in `validateBamRunPaths` in
`modules/helpers/path_validation.nf`; otherwise preflight will reject the inputs.

For pseudobulk, use `params.pseudobulk = true`, set `params.ncells`, and
arrange files like this for `params.ncells = 10000`:

```text
<pseudobulkDir>/
  my_dataset/
    sampleA/
      pseudobulk_n_barcodes_10000_1.bam
      pseudobulk_n_barcodes_10000_2.bam
```

The glob accepts other suffixes, but Conpair's pseudobulk pairing in
`main.nf` requires files ending in `_1.bam` and `_2.bam`; it compares the
two sets pairwise. Use the same naming convention for other donors.

For either BAM layout, check read groups and sample (`SM`) identifiers
before genotyping. The VCF sample identifiers produced by bcftools need
to be unique across inputs and include the corresponding modality (or
`pseudobulk`), since downstream VCF grouping uses those identifiers. The
pipeline renames BAM files but does not generally rename their read-group
sample identifiers before genotype calling.

### Existing VCF input: `main_from_genotype.nf`

Supply one coordinate-sorted, bgzip-compressed VCF per sample, with a
matching CSI index named `<filename>.vcf.gz.csi`. The pipeline merges the
files with `bcftools merge`, so sample identifiers must be unique across
files and records must be compatible with merging. The VCFs must contain
genotypes and depth information for the pipeline's `DP` filtering.

For real data, set `params.modalities` in `nextflow.config` to your two
modality names, as for the BAM workflow. For example:

```text
<vcfsDir>/
  my_dataset/
    real_data/ncells_null/
      sampleA_bulk.vcf.gz
      sampleA_bulk.vcf.gz.csi
      sampleA_single-cell.vcf.gz
      sampleA_single-cell.vcf.gz.csi
```

Each filename **and** the sample identifier inside its VCF should contain
the corresponding modality: individual files are grouped by filename,
while merged VCF samples are grouped by their internal identifiers. For
example, `bcftools query -l sampleA_bulk.vcf.gz` should list a unique
identifier such as `sampleA_bulk`. A shared donor identifier can still be
used to record the expected match across modalities.

For pseudobulk, use `params.pseudobulk = true`, set `params.ncells`, and
put files and their CSI indexes under:

```text
<vcfsDir>/my_dataset/pseudobulk/ncells_10000/
  sampleA_pseudobulk_1.vcf.gz
  sampleA_pseudobulk_1.vcf.gz.csi
  sampleA_pseudobulk_2.vcf.gz
  sampleA_pseudobulk_2.vcf.gz.csi
```

Use unique VCF sample identifiers containing `pseudobulk` for these files
too. `main.nf` publishes its individual genotypes under
`<outdir>/1a_individual_vcf/` in this VCF layout, but check their internal
sample names before using them as input to another run.

### Special VCF experiments

These workflows search `params.vcfsDir` using the layout above. Review their
workflow-specific parameters and `params.outdir` assignments in the entry
point before running. The missing-sample, double-sample, and uneven-size
workflows are pseudobulk-only and intentionally use the fixed
`["pseudobulk"]` modality label; `params.modalities` does not apply to them.
The pseudobulk-versus-single-cell workflow also has fixed modality labels:
it compares `pseudobulk` with `single-cell` and does not use
`params.modalities`. Set `params.pseudobulk = true` for the first three
experiments; for pseudobulk versus single-cell set it to `false` so that
the downstream tools use the real-data branch:

- Missing samples: `main_missing_samples.nf` uses `_1.vcf.gz` and
  `_2.vcf.gz`; `n_samples_to_remove` controls removal from set 2 and
  `repeat` supplies the random seed.
- Double samples: `main_double_samples.nf` uses sets `_1`, `_2` and `_3`;
  `n_double_samples` selects additional files from set 3 and `repeat`
  affects the seed.
- Uneven sizes: `main_uneven_pseudobulk_sizes.nf` selects `_1` from
  `ncells_${params.ncells}` and `_2` from `ncells_${params.ncells_2}`.
- Pseudobulk versus single-cell: `main_pseudobulk_vs_sc.nf` selects
  `_1.vcf.gz` pseudobulk files and `_single-cell.vcf.gz` real-data files.

Supply matching CSI indexes for all selected VCFs. These experiment
workflows run a subset of the tools and should be checked carefully with
your own data and sample naming.

## 4. Set read-depth thresholds

`params.read_depth` in `nextflow.config` contains entries by dataset and
modality. For an unknown dataset, the existing fallback has `default: 0`
and `upper: 1000`; this is not an automatically chosen threshold for your
data. Review depth distributions and add an entry with explicit thresholds
for both real-data modalities, for example:

```groovy
"my_dataset": [
    "bulk": 3,
    "single-cell": 2,
    "upper": 100,
    "default": 0
],
```

These numbers are illustrative, not recommendations. For pseudobulk the
helper uses the dataset's `default` threshold and `upper` cap. For real
data, set each modality explicitly: a missing modality currently falls back
to zero during filtering even if the dataset entry has a nonzero `default`.
Read-depth filtering uses VCF `DP`; `params.mpileup_max_depth` separately
limits depth during BAM genotype calling. Review the other tool-specific
settings in the [parameter notes](README.md#parameters).

## 5. Check assumptions and run

1. Confirm both modalities or pseudobulk sets are discovered using the
   patterns above; verify all files and indexes are readable on compute
   nodes. Check VCF identifiers with `bcftools query -l` and BAM read
   groups with `samtools view -H`.
2. The pipeline still uses sample names to assign samples to modalities, and
   applies the configured quality and read-depth filters.
3. Keep a separate mapping of expected donor matches if you want to assess
   benchmark accuracy. The pipeline does not accept a ground-truth sample
   sheet or automatically establish expected matches.
4. Choose an appropriately small initial dataset and run, for example:

   ```bash
   nextflow run main.nf -profile slurm -w nextflow -with-trace
   ```

   Substitute `main_from_genotype.nf` for VCF inputs, or `local` for the
   profile when appropriate. See the [launch instructions](README.md#running-the-pipeline)
   for Slurm setup and the optional `-resume` flag.
5. Check the discovered-file messages, `.nextflow.log`, trace report, task
   logs in the work directory, and per-tool directories in `params.outdir`.
   A completed workflow alone does not guarantee every tool succeeded:
   Vireo uses `errorStrategy 'ignore'`, so verify its expected outputs as
   well as the other tools' results.

#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { isValidRealDataModalities; validateBenchmarkParameters } from '../modules/helpers/modality_utils'
include {
    ensureDirectory; ensureFile; findPathMatches; requirePositiveInteger;
    throwValidationErrors; validateOutputDirectory
} from '../modules/helpers/path_validation'

include { PROCESS_METADATA } from './modules/process_metadata.nf'
include { PSEUDOBULK_BAM_TO_FASTQ } from './modules/pseudobulk_bam_to_fastq.nf'
include { NTSM_COUNT } from './modules/ntsm_count.nf'
include { NTSM_EVAL } from './modules/ntsm_eval.nf'

def validateNtsmPaths(params, modalities) {
    def errors = []
    ensureFile(params.sites, 'NTSM sites FASTA', errors)
    validateOutputDirectory(params.outdir, errors)

    if (params.pseudobulk) {
        def pattern = "${params.prepared_bam_dir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.bam"
        def bams = findPathMatches(pattern)
        if (bams.isEmpty()) {
            errors << "No pseudobulk BAM files matched ${pattern}"
        }
        bams.each { bam -> ensureFile(bam, 'Pseudobulk BAM', errors) }
    } else if (params.dataset == 'hgsoc-new') {
        def r1Pattern = "${params.pl_hgsoc}/*/*/HGSOC-*_R1_001.fastq.gz"
        def r2Pattern = "${params.pl_hgsoc}/*/*/HGSOC-*_R2_001.fastq.gz"
        def r1Files = findPathMatches(r1Pattern)
        def r2Files = findPathMatches(r2Pattern)
        def r1Paths = r1Files.collect { it.absolutePath } as Set
        def r2Paths = r2Files.collect { it.absolutePath } as Set

        if (r1Files.isEmpty()) errors << "No R1 FASTQ files matched ${r1Pattern}"
        if (r2Files.isEmpty()) errors << "No R2 FASTQ files matched ${r2Pattern}"
        (r1Files + r2Files).each { fastq ->
            ensureFile(fastq, 'FASTQ input', errors)
            if (!(fastq.parentFile.parentFile.name in modalities)) {
                errors << "FASTQ modality directory '${fastq.parentFile.parentFile.name}' is not in params.modalities: ${fastq}"
            }
        }
        modalities.each { mod ->
            if (!r1Files.any { it.parentFile.parentFile.name == mod }) {
                errors << "No R1 FASTQ files found for configured modality '${mod}'."
            }
        }
        r1Files.each { fastq ->
            def mate = fastq.absolutePath.replaceFirst(/_R1_001\.fastq\.gz$/, '_R2_001.fastq.gz')
            if (!r2Paths.contains(mate)) errors << "Missing R2 FASTQ mate for ${fastq}"
        }
        r2Files.each { fastq ->
            def mate = fastq.absolutePath.replaceFirst(/_R2_001\.fastq\.gz$/, '_R1_001.fastq.gz')
            if (!r1Paths.contains(mate)) errors << "Missing R1 FASTQ mate for ${fastq}"
        }
    } else {
        ensureFile("${params.projectDir}/process_metadata.py", 'NTSM metadata helper script', errors)
        modalities.each { mod ->
            ensureFile(
                "${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${mod}_${params.dataset}_${params.read_type}.csv",
                "Metadata for modality '${mod}'",
                errors
            )
            ensureDirectory("${params.petaLibrary}/raw_reads/${params.dataset}/${mod}", "FASTQ directory for modality '${mod}'", errors)
        }
    }

    throwValidationErrors(errors)
}

def validateMappedFastqs(row, modality, readType) {
    def errors = []
    def sample = row.'Sample Name'?.toString()
    if (!sample || !(sample ==~ /[A-Za-z0-9._-]+/) || sample in ['.', '..']) {
        errors << "Metadata for modality '${modality}' has an invalid sample name: '${sample}'"
    }
    def runKeys = row.keySet().findAll { it ==~ /R[12]_path_run_\d+/ }
        .collect { it.replaceFirst(/^R[12]_path_/, '') }.unique()
    def runsWithReads = 0
    runKeys.each { run ->
        def r1 = row["R1_path_${run}"]?.toString()?.trim()
        def r2 = row["R2_path_${run}"]?.toString()?.trim()
        if (r1 || r2) {
            runsWithReads++
            if (!r1) errors << "Missing R1 path for sample '${sample}', ${run}."
            else ensureFile(r1, "R1 FASTQ for sample '${sample}'", errors)
            if (readType == 'paired') {
                if (!r2) errors << "Missing R2 path for sample '${sample}', ${run}."
                else ensureFile(r2, "R2 FASTQ for sample '${sample}'", errors)
            }
        }
    }
    if (!runsWithReads) errors << "No FASTQ paths found for sample '${sample}' in modality '${modality}'."
    return errors
}

// Collect all columns matching a run prefix (e.g., R1_path_run_0, R1_path_run_1)
// and return them in numeric run order.
def collectRunsDynamic(row, prefix) {
    // collect keys and build a plain-string regex (no slashy regex interpolation)
    def keys   = row.keySet().toList()
    def regex  = "${prefix}_run_\\d+"        // note: double-escape for \d inside a GString
    def matched = keys.findAll { k -> k ==~ regex }
    def sorted = matched.sort { a, b ->
        (a =~ /\d+$/)[0].toInteger() <=> (b =~ /\d+$/)[0].toInteger()
    }
    def values = sorted.collect { k -> row[k]?.toString()?.trim() }.findAll { value -> value }
    def files = values.collect { value -> file(value) }

    return files
}

// Single-end runs have no real R2 files; lazily create an empty gzipped
// placeholder so NTSM_COUNT's required R2 path input is always satisfied.
def emptyR2Placeholder() {
    def placeholder = file("${workDir}/ntsm_empty_R2.fastq.gz")
    if (!placeholder.exists()) {
        placeholder.parent.mkdirs()
        new java.util.zip.GZIPOutputStream(new FileOutputStream(placeholder.toString())).close()
    }
    return placeholder
}

workflow {
    validateBenchmarkParameters(params.dataset, params.pseudobulk, params.ncells)
    requirePositiveInteger(params.threads, 'params.threads')
    if (!(params.read_type in ['single', 'paired'])) {
        throw new IllegalArgumentException("params.read_type must be 'single' or 'paired'.")
    }
    if (!params.pseudobulk && (!(params.datatype instanceof String) ||
        !(params.datatype ==~ /[A-Za-z0-9._-]+/) || params.datatype in ['.', '..'])) {
        throw new IllegalArgumentException("params.datatype must be a safe directory name using only letters, numbers, '.', '_' or '-'.")
    }
    if (!params.pseudobulk && !isValidRealDataModalities(params.modalities)) {
        throw new IllegalArgumentException(
            "Real-data runs require exactly two distinct, safe modality names in params.modalities; neither may contain the other."
        )
    }
    modalities = params.pseudobulk ? ['pseudobulk'] : params.modalities
    validateNtsmPaths(params, modalities)

    if (params.pseudobulk) {
        // For pseudobulk experiments, we need to convert pseudobulk bams to fastqs
        pseudobulk_bam_pattern = "${params.prepared_bam_dir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.bam"
        pseudobulk_bams = channel
            .fromPath(pseudobulk_bam_pattern)
            .map { bam -> tuple(bam.baseName.toString(), bam) }
            .view { x -> "Pseudobulk BAM files: ${x}" }

        pseudobulk_fastqs = PSEUDOBULK_BAM_TO_FASTQ(pseudobulk_bams)

        sample_fastqs = pseudobulk_fastqs
            .map { sample_name, r1_fastq, r2_fastq ->
                tuple('pseudobulk', sample_name, [r1_fastq], [r2_fastq])
            }
            .view { row -> "Pseudobulk sample fastqs: ${row}" }
    }
    else if (params.dataset == "hgsoc-new") {
        // For the new HGSOC, we do not need to process metadata CSVs, as the fastq files are already organized in a simple directory structure.
        fastq_r1_samples = channel
            .fromPath("${params.pl_hgsoc}/*/*/HGSOC-*_R1_001.fastq.gz")
            .map { r1_fastq ->
                tuple(r1_fastq.parent.parent.name, r1_fastq.parent.name.take(10), r1_fastq)
            }
            .groupTuple(by: [0, 1], sort: true)

        fastq_r2_samples = channel
            .fromPath("${params.pl_hgsoc}/*/*/HGSOC-*_R2_001.fastq.gz")
            .map { r2_fastq ->
                tuple(r2_fastq.parent.parent.name, r2_fastq.parent.name.take(10), r2_fastq)
            }
            .groupTuple(by: [0, 1], sort: true)

        sample_fastqs = fastq_r1_samples
            .join(fastq_r2_samples, by: [0, 1])
    }
    else {
        // For other datasets, we need to process metadata and fastq files.
        // Key metadata and fastq dir channels by modality, then join so each
        // metadata file is always paired with its own modality's fastq dir.
        metadata_ch = channel.fromList(modalities)
            .map { mod ->
                tuple(mod, file("${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${mod}_${params.dataset}_${params.read_type}.csv"))
            }
        fastq_dir_ch = channel.fromList(modalities)
            .map { mod ->
                tuple(mod, file("${params.petaLibrary}/raw_reads/${params.dataset}/${mod}", type: 'dir').toAbsolutePath())
            }
        metadata_ch.view { x -> "Metadata files: ${x}"}
        fastq_dir_ch.view { x -> "Fastq files: ${x}" }

        // join on modality so each metadata file is paired with its own fastq dir
        // in a single 3-element tuple, matching PROCESS_METADATA's input declaration
        process_metadata_ch = metadata_ch.join(fastq_dir_ch)

        mapping_csv = PROCESS_METADATA(process_metadata_ch)

        // elem: 1 tells splitCsv the CSV path lives at tuple index 1, keeping modality at index 0
        sample_fastqs = mapping_csv.map { mod, mapping ->
                def firstRow = mapping.toFile().withReader('UTF-8') { reader ->
                    reader.readLine()
                    reader.readLine()
                }
                if (!firstRow?.trim()) {
                    throw new IllegalStateException("Metadata for modality '${mod}' produced no sample mappings: ${mapping}")
                }
                tuple(mod, mapping)
            }
            .splitCsv(header: true, elem: 1)
                .map { mod, row ->
                    def errors = validateMappedFastqs(row, mod, params.read_type)
                    if (errors) {
                        throw new IllegalStateException("FASTQ validation failed:\n - ${errors.unique().join('\n - ')}")
                    }
                    def r2_paths = params.read_type == 'single' ? [emptyR2Placeholder()] : collectRunsDynamic(row, 'R2_path')
                    tuple(
                        mod,
                        row.'Sample Name',
                        collectRunsDynamic(row, 'R1_path'),
                        r2_paths
                    )
                }
            .view { row -> "Sample fastq dirs: ${row}" }
    } 
    
    sample_counts = NTSM_COUNT(sample_fastqs)
    NTSM_EVAL(sample_counts.collect(), modalities.collect())
}

def ensureDirectory(path, label, errors) {
    if (path == null || path.toString().trim().isEmpty()) {
        errors << "${label} must be set to a directory path."
        return
    }
    def directory = new File(path?.toString() ?: '')
    if (!directory.isDirectory() || !directory.canRead()) {
        errors << "${label} must be an existing, readable directory: ${directory}"
    }
}

def ensureFile(path, label, errors) {
    def file = new File(path?.toString() ?: '')
    if (!file.isFile() || !file.canRead() || file.length() == 0) {
        errors << "${label} must be an existing, non-empty, readable file: ${file}"
    }
}

def ensureExecutable(path, label, errors) {
    def file = new File(path?.toString() ?: '')
    if (!file.isFile() || !file.canRead() || !file.canExecute()) {
        errors << "${label} must be an existing, readable, executable file: ${file}"
    }
}

def ensureAnyFile(paths, label, errors) {
    def files = paths.collect { new File(it?.toString() ?: '') }
    if (!files.any { it.isFile() && it.canRead() && it.length() > 0 }) {
        errors << "${label} must exist and be readable; checked: ${files.join(', ')}"
    }
}

def findPathMatches(pattern) {
    def absolutePattern = new File(pattern.toString()).toPath().toAbsolutePath().normalize().toString()
    def wildcardPositions = ['*', '?', '{', '[']
        .collect { absolutePattern.indexOf(it) }
        .findAll { it >= 0 }

    if (wildcardPositions.isEmpty()) {
        def exactFile = new File(absolutePattern)
        return exactFile.exists() ? [exactFile] : []
    }

    def firstWildcard = wildcardPositions.min()
    def prefix = absolutePattern.substring(0, firstWildcard)
    def rootEnd = prefix.lastIndexOf(File.separator)
    def root = new File(rootEnd > 0 ? prefix.substring(0, rootEnd) : File.separator)
    if (!root.isDirectory()) {
        return []
    }

    def matcher = java.nio.file.FileSystems.getDefault().getPathMatcher("glob:${absolutePattern}")
    def matches = []
    def paths
    try {
        paths = java.nio.file.Files.walk(
            root.toPath(),
            java.nio.file.FileVisitOption.FOLLOW_LINKS
        )
        paths.forEach { path ->
            def candidate = path.toAbsolutePath().normalize()
            if (matcher.matches(candidate) && java.nio.file.Files.isRegularFile(candidate)) {
                matches << candidate.toFile()
            }
        }
    } catch (java.io.IOException e) {
        throw new IllegalStateException("Could not search for input files matching ${pattern}: ${e.message}", e)
    } finally {
        paths?.close()
    }
    matches.sort { a, b -> a.path <=> b.path }
}

def validateVcfFiles(files, label, errors) {
    files.each { vcf ->
        ensureFile(vcf, "${label} VCF", errors)
        ensureFile("${vcf}.csi", "${label} CSI index for ${vcf.name}", errors)
    }
}

def validateVcfRunPaths(params, tools, inputPatterns, modalities = []) {
    def errors = []
    ensureDirectory(params.projectDir, 'params.projectDir', errors)
    validateOutputDirectory(params.outdir, errors)
    validateToolPaths(params, tools, errors)

    def inputs = [:]
    inputPatterns.each { label, pattern ->
        def patterns = pattern instanceof List ? pattern : [pattern]
        def files = patterns.collectMany { findPathMatches(it) }.unique { it.path }
        if (files.isEmpty()) {
            errors << "${label}: no input files matched ${patterns.join(' or ')}"
        }
        validateVcfFiles(files, label, errors)
        inputs[label] = files
    }

    modalities.each { modality ->
        def files = inputs.values().flatten()
        if (!files.any { it.name.contains(modality) }) {
            errors << "No VCF filename contains configured modality '${modality}'. Check params.modalities and the input filenames."
        }
    }

    throwValidationErrors(errors)
    inputs
}

def validateBamRunPaths(params, tools) {
    def errors = []
    ensureDirectory(params.projectDir, 'params.projectDir', errors)
    validateOutputDirectory(params.outdir, errors)
    validateToolPaths(params, tools, errors)

    def bamGroups = [:]
    if (params.pseudobulk) {
        def pattern = "${params.pseudobulkDir}/${params.dataset}/*/pseudobulk_n_barcodes_${params.ncells}_*.bam"
        bamGroups.pseudobulk = findPathMatches(pattern)
        if (bamGroups.pseudobulk.isEmpty()) {
            errors << "No pseudobulk BAM files matched ${pattern}"
        }
        if (!bamGroups.pseudobulk.any { it.name ==~ /.*_1\.bam$/ } ||
            !bamGroups.pseudobulk.any { it.name ==~ /.*_2\.bam$/ }) {
            errors << "Pseudobulk BAM inputs must include both *_1.bam and *_2.bam files for Conpair pairing."
        }
    } else {
        params.modalities.each { modality ->
            def pattern = "${params.bamsDir}/${params.dataset}/${modality}/*/Aligned.sortedByCoord.out.bam"
            bamGroups[modality] = findPathMatches(pattern)
            if (bamGroups[modality].isEmpty()) {
                errors << "No BAM files for modality '${modality}' matched ${pattern}"
            }
        }
    }
    bamGroups.values().flatten().each { bam ->
        ensureFile(bam, 'BAM input', errors)
    }

    def fasta = "${params.refGenome}/fasta/genome.fa"
    ensureFile(fasta, 'Reference FASTA', errors)
    ensureFile("${fasta}.fai", 'Reference FASTA index', errors)
    ensureAnyFile(
        ["${params.refGenome}/fasta/genome.dict", "${fasta}.dict"],
        'Reference sequence dictionary',
        errors
    )
    ensureFile(params.snp_vcf, 'SNP reference VCF', errors)

    throwValidationErrors(errors)
    bamGroups
}

def requirePositiveInteger(value, label) {
    def numericValue = value instanceof Number ? value.doubleValue() : Double.NaN
    if (!Double.isFinite(numericValue) || numericValue <= 0 || numericValue != Math.floor(numericValue)) {
        throw new IllegalArgumentException("${label} must be a positive integer.")
    }
}

def validateOutputDirectory(path, errors) {
    if (path == null || path.toString().trim().isEmpty()) {
        errors << "params.outdir must be set to a directory path."
        return
    }
    def directory = new File(path?.toString() ?: '')
    if (directory.exists() && !directory.isDirectory()) {
        errors << "params.outdir must be a directory path: ${directory}"
        return
    }
    if (!directory.exists() && !directory.mkdirs()) {
        errors << "Could not create params.outdir: ${directory}"
        return
    }
    try {
        def probe = java.nio.file.Files.createTempFile(directory.toPath(), '.nextflow-preflight-', '.tmp')
        java.nio.file.Files.delete(probe)
    } catch (java.io.IOException e) {
        errors << "params.outdir is not writable: ${directory} (${e.message})"
    }
}

def validateToolPaths(params, tools, errors) {
    def project = params.projectDir
    if (tools.contains('vireo')) {
        ensureFile("${project}/modules/vireo.py", 'Vireo helper script', errors)
    }
    if (tools.contains('fingerprints')) {
        ensureFile("${project}/modules/hg38_chr.map", 'CrosscheckFingerprints haplotype map', errors)
    }
    if (tools.contains('hysys')) {
        ensureExecutable("${params.HYSYS}/HaveYouSwappedYourSamples.sh", 'HYSYS launcher', errors)
        ensureFile("${params.HYSYS}/concordance_script.py", 'HYSYS concordance script', errors)
        ensureFile("${params.HYSYS}/model_analysis.py", 'HYSYS model script', errors)
    }
    if (tools.contains('ngscheckmate')) {
        ensureFile("${params.NGS_CHECKMATE}/ncm.py", 'NGSCheckMate launcher', errors)
        ensureFile("${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed", 'NGSCheckMate SNP BED', errors)
    }
    if (tools.contains('timeattackgencomp')) {
        ensureFile("${params.TIMEATTACKGENCOMP}/compare_simple.pl", 'TimeAttackGenComp comparison script', errors)
        ensureFile(params.timeattackgencomp_target_bed, 'TimeAttackGenComp target BED', errors)
    }
    if (tools.contains('conpair')) {
        ensureFile("${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py", 'Conpair pileup script', errors)
        ensureFile("${params.CONPAIR}/scripts/verify_concordance.py", 'Conpair verification script', errors)
        ensureDirectory("${params.CONPAIR}/modules", 'Conpair modules directory', errors)
        ensureFile(
            "${project}/modules/helpers/conpair_normalize_pileup.py",
            'Conpair pileup normalizer',
            errors
        )
        def defaultMarker = "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover"
        ensureFile(params.conpair_marker_bed ?: "${defaultMarker}.bed", 'Conpair marker BED', errors)
        ensureFile(params.conpair_marker_txt ?: "${defaultMarker}.txt", 'Conpair marker TXT', errors)
    }
}

def throwValidationErrors(errors) {
    if (!errors.isEmpty()) {
        throw new IllegalStateException("Preflight validation failed:\n - ${errors.unique().join('\n - ')}")
    }
}

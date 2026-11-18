def isValidRealDataModalities(modalities) {
    modalities instanceof List &&
        modalities.size() == 2 &&
        modalities.every { it instanceof String && it ==~ /[A-Za-z0-9._-]+/ && !(it in ['.', '..']) } &&
        modalities.toSet().size() == 2 &&
        !modalities[0].contains(modalities[1]) &&
        !modalities[1].contains(modalities[0])
}

def validateBenchmarkParameters(dataset, pseudobulk, ncells) {
    if (!(dataset instanceof String) || !(dataset ==~ /[A-Za-z0-9._-]+/) || dataset in ['.', '..']) {
        throw new IllegalArgumentException("params.dataset must be a safe directory name using only letters, numbers, '.', '_' or '-'.")
    }

    if (!(pseudobulk instanceof Boolean)) {
        throw new IllegalArgumentException("params.pseudobulk must be a boolean (true or false).")
    }

    if (pseudobulk) {
        def cellCount = ncells instanceof Number ? ncells.doubleValue() : Double.NaN
        if (!Double.isFinite(cellCount) || cellCount <= 0 || cellCount != Math.floor(cellCount)) {
            throw new IllegalArgumentException("params.ncells must be a positive integer when params.pseudobulk is true.")
        }
    }
}

#!/usr/bin/env nextflow

include { buildModalityOutputContext } from './helpers/read_depth_utils'

process CONPAIR_PILEUP {
    conda "${params.conda}/conpair"
    errorStrategy 'ignore'

    input:
        tuple val(bam_key), val(pileup_id), path(bam)

    output:
        tuple val(bam_key), path("${pileup_id}_pileup.txt"), emit: pileup

    script:
        def conpairMarkerBed = params.conpair_marker_bed ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
        def conpairMarkerTxt = params.conpair_marker_txt ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"

        """
        set -euo pipefail

        export CONPAIR_DIR="${params.CONPAIR}"
        export GATK_JAR="${params.conda}/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar"
        export PYTHONPATH="${params.CONPAIR}/modules"

        conpair_marker_bed="${conpairMarkerBed}"
        conpair_marker_txt="${conpairMarkerTxt}"
        conpair_normalize_pileup_py="${params.projectDir}/modules/helpers/conpair_normalize_pileup.py"
        pileup="${pileup_id}_pileup.txt"
        pileup_raw="${pileup_id}_pileup.raw.txt"

        if [ ! -s "\$conpair_marker_bed" ]
        then
            echo "ERROR: Conpair BED marker file not found: \$conpair_marker_bed" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_marker_txt" ]
        then
            echo "ERROR: Conpair TXT marker file not found: \$conpair_marker_txt" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_normalize_pileup_py" ]
        then
            echo "ERROR: Conpair pileup normalizer not found: \$conpair_normalize_pileup_py" >&2
            exit 1
        fi

        fasta_contig="\$(grep -m1 '^>' "${params.refGenome}/fasta/genome.fa" | sed 's/^>//' | awk '{print \$1}')"
        marker_contig="\$(awk 'NF {print \$1; exit}' "\$conpair_marker_txt")"
        remove_chr_opt=""

        if [[ "\$fasta_contig" == chr* && "\$marker_contig" != chr* ]]
        then
            remove_chr_opt="--remove_chr_prefix"
        elif [[ "\$fasta_contig" != chr* && "\$marker_contig" == chr* ]]
        then
            echo "ERROR: Reference FASTA contigs do not use 'chr' prefix, but selected Conpair markers do: \$conpair_marker_txt" >&2
            echo "Set --conpair_marker_bed/--conpair_marker_txt to marker files that match your reference naming." >&2
            exit 1
        fi

        python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \\
                -B ${bam} \\
                -R ${params.refGenome}/fasta/genome.fa \\
                -M \$conpair_marker_bed \\
                \$remove_chr_opt \\
                -O "\$pileup_raw"

        python "\$conpair_normalize_pileup_py" "\$pileup_raw" "\$pileup.tmp"

        python - "\$pileup.tmp" "\$pileup" <<'PY'
import sys

in_path, out_path = sys.argv[1], sys.argv[2]

with open(in_path, "r", encoding="utf-8", errors="replace") as fin, open(out_path, "w", encoding="utf-8") as fout:
    for line in fin:
        if not line.strip() or line.startswith("[REDUCE RESULT]"):
            continue

        fields = line.rstrip("\\n").split("\\t")
        if len(fields) < 5:
            fields = line.rstrip("\\n").split()
        if len(fields) < 5:
            continue

        chrom, pos, ref = fields[0], fields[1], fields[2].upper()
        bases = fields[3].upper()
        quals = fields[4]
        filtered_pairs = [(base, qual) for base, qual in zip(bases, quals) if base in {"A", "C", "G", "T"}]
        if not filtered_pairs:
            continue

        out_bases = "".join(base for base, _ in filtered_pairs)
        out_quals = "".join(qual for _, qual in filtered_pairs)
        if len(out_bases) != len(out_quals):
            continue

        fout.write(f"{chrom}\\t{pos}\\t{ref}\\t{out_bases}\\t{out_quals}\\n")
PY

        rm -f "\$pileup.tmp"

        if [ ! -s "\$pileup" ]
        then
            echo "ERROR: Sanitized Conpair pileup is empty: \$pileup" >&2
            exit 1
        fi
        """
}

process CONPAIR_VERIFY {
    conda "${params.conda}/conpair"
    publishDir "${params.outdir}/2b_conpair", mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple val(sample_1), path(pileup_1), val(sample_2), path(pileup_2), val(mod_1), val(mod_2)

    output:
        path("${params.dataset}/**/${sample_1}_vs_${sample_2}_concordance.txt")

    script:
        def conpairMarkerTxt = params.conpair_marker_txt ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
        def conpairMinCov = params.conpair_min_cov ?: 3
        def conpairMinMapQual = params.conpair_min_map_qual ?: 0
        def conpairMinBaseQual = params.conpair_min_base_qual ?: 10
        def outputContext = buildModalityOutputContext(
            params.dataset,
            params.pseudobulk,
            params.ncells,
            [mod_1, mod_2]
        )

        """
        set -euo pipefail

        export CONPAIR_DIR="${params.CONPAIR}"
        export GATK_JAR="${params.conda}/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar"
        export PYTHONPATH="${params.CONPAIR}/modules"

        conpair_marker_txt="${conpairMarkerTxt}"
        conpair_min_cov=${conpairMinCov}
        conpair_min_map_qual=${conpairMinMapQual}
        conpair_min_base_qual=${conpairMinBaseQual}
        output_location="${outputContext.experimentPath}"
        output_file="\${output_location}/${sample_1}_vs_${sample_2}_concordance.txt"

        if [ ! -s "\$conpair_marker_txt" ]
        then
            echo "ERROR: Conpair TXT marker file not found: \$conpair_marker_txt" >&2
            exit 1
        fi

        pileup_contig_1="\$(awk 'NF {print \$1; exit}' ${pileup_1})"
        pileup_contig_2="\$(awk 'NF {print \$1; exit}' ${pileup_2})"
        if [[ -z "\$pileup_contig_1" || -z "\$pileup_contig_2" ]]
        then
            echo "ERROR: Empty pileup detected for one or both samples: ${pileup_1}, ${pileup_2}" >&2
            exit 1
        fi

        mkdir -p "\$(dirname "\${output_file}")"

        python ${params.CONPAIR}/scripts/verify_concordance.py \\
            -T ${pileup_1} \\
            -N ${pileup_2} \\
            -M \$conpair_marker_txt \\
            -C \$conpair_min_cov \\
            -Q \$conpair_min_map_qual \\
            -B \$conpair_min_base_qual \\
            -O "\${output_file}"

        if [ ! -s "\${output_file}" ]
        then
            echo "WARN: No concordance output for ${sample_1} vs ${sample_2} at C=\$conpair_min_cov,Q=\$conpair_min_map_qual,B=\$conpair_min_base_qual; retrying with C=1,Q=0,B=0" >&2
            python ${params.CONPAIR}/scripts/verify_concordance.py \\
                -T ${pileup_1} \\
                -N ${pileup_2} \\
                -M \$conpair_marker_txt \\
                -C 1 \\
                -Q 0 \\
                -B 0 \\
                -O "\${output_file}"
        fi

        if [ ! -s "\${output_file}" ]
        then
            echo "ERROR: Conpair did not produce concordance output for ${sample_1} vs ${sample_2} even after fallback thresholds." >&2
            echo "Check pileup depth and marker/reference compatibility." >&2
            exit 1
        fi
        """
}
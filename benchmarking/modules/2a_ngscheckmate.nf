#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process NGSCHECKMATE {
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/2a_ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf)
        val(modalities)
    
    output:
        path("${params.dataset}/**/output_all.txt")
        path("${params.dataset}/**/output_matched.txt")
        path("${params.dataset}/**/output_output_corr_matrix.txt")

    script:
        def readDepthContext = buildReadDepthContext(
            params.read_depth,
            params.dataset,
            params.pseudobulk,
            modalities,
            params.ncells
        )

        def outputPath = "${readDepthContext.experimentPath}/read_depth_${params.pseudobulk ? readDepthContext.pseudobulkReadDepth : readDepthContext.modalityReadDepthTag}"

        """
        set -euo pipefail

        output_location="${outputPath}"
        mkdir -p \$output_location

        # Unzip and list
        for file in ${vcf}
        do  
            base_name=\$(basename "\$file")
            base_name=\${base_name%.gz}

            gunzip -c \${file} > \${output_location}/\${base_name}
                
            echo "\${output_location}/\${base_name}" >> \${output_location}/vcf_list.txt
        done
        
        # Run NGSCheckMate in VCF mode
        python ${params.NGS_CHECKMATE}/ncm.py \\
            -V \\
            -l \${output_location}/vcf_list.txt \\
            -bed ${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed \\
            -O \${output_location}
        """
}
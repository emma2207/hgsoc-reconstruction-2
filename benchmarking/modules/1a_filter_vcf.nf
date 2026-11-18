#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process FILTER_VCF {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_vcf", mode: 'copy'
    
    input:
        path(vcf_file)
        val(modalities)

    output:
        path("${params.dataset}/**/*_individual_variants.vcf.gz"), emit: individual_vcfs
        path("${params.dataset}/**/*_modality_variants.vcf.gz")
        path("${params.dataset}/**/filtered_variants_*_rd_*.vcf.gz"), emit: modality_vcfs

    script:
        def readDepthContext = buildReadDepthContext(
            params.read_depth,
            params.dataset,
            params.pseudobulk,
            modalities
        )

        def modalityReadDepthTag = readDepthContext.modalityReadDepthTag
        def pseudobulkReadDepth = readDepthContext.pseudobulkReadDepth
        def readDepthPairs = readDepthContext.readDepthPairs
        def realDataNcellsPath = readDepthContext.realDataNcellsPath

        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/${realDataNcellsPath}/read_depth_${modalityReadDepthTag}" 
        fi
        mkdir -p \$output_location

        READ_DEPTH_MAP='${readDepthPairs}'
        declare -A RD_BY_MOD
        IFS=',' read -ra RD_PAIRS <<< "\$READ_DEPTH_MAP"
        for pair in "\${RD_PAIRS[@]}"
        do
            key="\${pair%%:*}"
            value="\${pair#*:}"
            RD_BY_MOD["\$key"]="\$value"
        done
        DEFAULT_RD="\${RD_BY_MOD[default]:-0}"

        # List samples by modality
        for sample in \$(bcftools query -l $vcf_file)
        do  
            for mod in ${modalities.join(' ')}
            do  
                if [[ "\$sample" ==  *"\$mod"* ]]
                then
                    echo "\$sample"
                    echo "\$sample" >> "\${output_location}/\${mod}_samples_vcf.txt"
                fi
            done
        done

        echo "Made lists"

        # Save variants separately for different modalities
        for mod in ${modalities.join(' ')}
        do
            echo "\$mod"
            mod_rd="\${RD_BY_MOD[\$mod]:-\$DEFAULT_RD}"

            bcftools view -S "\${output_location}/\${mod}_samples_vcf.txt" \
                -Oz -o "\${output_location}/\${mod}_modality_variants.vcf.gz" \
                $vcf_file
            echo "split"
            bcftools index \${output_location}/\${mod}_modality_variants.vcf.gz

            nr_samples=\$( bcftools query -l \${output_location}/\${mod}_modality_variants.vcf.gz | wc -l )
            read_depth_cutoff=\$(( mod_rd * nr_samples ))

            all_variants_output="\${output_location}/filtered_variants_\${mod}_rd_\${mod_rd}.vcf.gz"

            # Filter variants
            bcftools view \
                -Oz \
                -i "QUAL>=20 && DP>=\${read_depth_cutoff}" \
                -o \${all_variants_output} \
                \${output_location}/\${mod}_modality_variants.vcf.gz

            # Index the filtered VCF
            bcftools index \${all_variants_output}
        done
        
        echo "Finished splitting by modalities!"

        # Save variants individually
        for sample in \$(bcftools query -l $vcf_file)
        do  
            echo "\$sample"
            sample_id="\${sample%????}"
            echo \$sample_id

            sample_mod="default"
            for mod in ${modalities.join(' ')}
            do
                if [[ "\$sample" == *"\$mod"* ]]
                then
                    sample_mod="\$mod"
                    break
                fi
            done
            sample_rd="\${RD_BY_MOD[\$sample_mod]:-\$DEFAULT_RD}"

            bcftools view \
                -c1 -Oz -s \$sample \
                -i "QUAL>=20 && DP>=\${sample_rd}" \
                -o "\${output_location}/\${sample_id}_individual_variants.vcf.gz" \
                $vcf_file

            bcftools index \${output_location}/\${sample_id}_individual_variants.vcf.gz
        done

        echo "Finished splitting output by sample"
        """
}
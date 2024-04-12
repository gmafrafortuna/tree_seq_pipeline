if len(allFiles) != 1:
    rule merge:
        input:
            vcf = [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)],
            idx = [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz.csi', file=allFiles)],
        output:
            vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz',
        params:
            files_path = f'{vcfdir}/{{chromosome}}/'
        conda: "bcftools"
        threads: 32
        resources: cpus=32, mem_mb=2048000, time_min=1200
        benchmark: 'benchmarks/{chromosome}.merge.benchmark.txt'
        shell:
            r"""
            printf '%s\n' {params.files_path}*.gz > files.merge
            
            bcftools merge \
                -m snps \
                -l files.merge \
                --threads {threads} -O z -o {output.vcf} 
            """

    rule index_merged:
        input: f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz',
        output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz.csi'
        conda: 'bcftools'
        shell:
            r"""
            bcftools index -f {input}
            """

rule remove_missing_indels:
    input: 
        vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz',
        index = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz.csi'
    output:    
        vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.recode.vcf.gz',
        index = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.recode.vcf.gz.csi'
    params:
        prefix = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged'
    conda: 'bcftools'
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
            --remove-indels \
            --max-missing 0.1 \
            --remove-filtered-all \
            --recode \
            --out {params.prefix}

        bgzip {params.prefix}.recode.vcf
        bcftools index -f {output.vcf}  

        rm  {params.prefix}.recode.vcf 
        """   

# rule summary:
#     input: 
#         vcf_file = rules.filter_merged_vcf.output.vcf,
#         index = rules.filter_merged_vcf.output.index,
#     output:
#         summary = f'{vcfdir}/{{chromosome}}/{{chromosome}}_summary.vchk'
#     params:
#         plots = f'{vcfdir}/{{chromosome}}/{{chromosome}}_summary/'
#     conda: 'bcftools'
#     benchmark: 'benchmarks/{chromosome}.summary.benchmark.txt'
#     shell:
#         r"""
#         bcftools stats {input.vcf_file} > {output.summary}
#         plot-vcfstats -p {params.plots} {output.summary}
#         """
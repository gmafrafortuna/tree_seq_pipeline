file_to_phase = rules.remove_missing_indels.output.vcf
file_to_phase_idx = rules.remove_missing_indels.output.index

phased_output = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz'
phased_output_idx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz.csi'

if config['ploidy'] == 1:
    rule rename_phased:
        input:
            vcf = file_to_phase,
            idx = file_to_phase_idx,
        output:
            vcf = temp(phased_output),
            idx = temp(phased_output_idx)
        log: 'logs/Rename_phased_{chromosome}.log'
        resources: cpus=1, mem_mb=32000, time_min=60
        shell:
            """
            if [ -h {input.vcf} ]; then
                ln -s $( realpath {input.vcf} ) {output.vcf}
                ln -s $( realpath {input.vcf} ).csi {output.vcf}
            else
                ln -s {input.vcf} {output.vcf}
                ln -s {input.idx} {output.idx}
            fi
            """
else:
    rule phase:
        input:
            vcf = file_to_phase,
            idx = file_to_phase_idx,
        output:
            vcf = phased_output,
        params:
            map = f'{mapdir}/{{chromosome}}.txt',
        threads: 32
        resources: cpus = 32, mem_mb = 2048, time_min = 7200
        conda: 'shapeit4am'
        benchmark:
            'benchmarks/{chromosome}.shapeit.benchmark.txt'
        shell:
            r"""
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            shapeit4 --input {input.vcf} \
                                --map {params.map} \
                                --region $chr \
                                --output {output.vcf} \
                                --thread {threads} \
                                --pbwt-depth 8 \
                                --sequencing \
                                --log {wildcards.chromosome}_phased.log                
            """         

    rule index_phased: 
        input:
            rules.phase.output.vcf
        output:
            idx = phased_output_idx
        conda: 'bcftools'
        benchmark:
            'benchmarks/{chromosome}.bcftools_index.bechmark.txt'
        shell:
            r"""
            bcftools index -f {input}
            """
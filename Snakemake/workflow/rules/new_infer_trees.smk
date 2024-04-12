# Convert vcf to samples file
# generate ancestos samples based only on sites with ancestral allele
# information
# truncate ancestors to remove those spanning too many sites
# generate ancestors trees
# match ancestors and samples to infer trees
# final output is a simplified tree sequence

rule prepare_sample_file:
    input:
        vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf.gz',
        meta = config['meta']
    output: f'{projdir}/Tsinfer/samples/{{chromosome}}.samples'
    params:
        chrLength = lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])]
    conda: "HLab_tsinfer"
    threads: 32
    resources: cpus = 32, mem_mb = 768000, time_min = 200
    benchmark:
            'benchmarks/{chromosome}_new.prep_sample.benchmark.txt'
    shell:
        r"""
        ulimit -v
        
        python scripts/generate_samples.py \
                {input.vcf} {input.meta} {output} {params.chrLength}
        
        """

rule generate_ancestors:
    input: 
        samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.samples',
        sites_to_exclude = f'{vcfdir}/{{chromosome}}/{{chromosome}}.exclude',
    output:
        ancestors_samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors',
    conda: 'HLab_tsinfer'
    threads: 64
    resources: cpus = 64, mem_mb = 2048000, time_min = 1680
    benchmark:
            'benchmarks/{chromosome}_generate_ancestors.benchmark.txt'
    shell:
        r"""
        # Check amount of memory (in kbytes) as seen by the job
        ulimit -v
        python scripts/generate_ancestors.py {input.samples} {input.sites_to_exclude} {output.ancestors_samples} {threads}
        """
# tsinfer generate-ancestors {input} \
#                 --num-threads {threads} \
#                 --num-flush-threads {threads} \
#                 --progress \
#                 --verbosity \
#                 --log-section tsinfer.inference


rule truncate_ancestors:
    input:
        ancestors_samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors',
    output: 
        truncated_ancestors = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors.truncated',
   # params:
   #     upper = 0.4,
   #     lower = 0.2,
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark:
            'benchmarks/{chromosome}_truncate.benchmark.txt'
    shell:
        r"""
        ulimit -v
        python scripts/truncate.py {input} {output}
        """

rule ancestors_trees:
    input:
        samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.samples',
        truncated_ancestors = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors.truncated',
    output: 
        ancestors_trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.atrees',
    params:
        re = '1.1e-8',
        mima = '1'
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus = 32, mem_mb = 2048000, time_min = 1680
    benchmark:
            'benchmarks/{chromosome}_ancestors_trees.benchmark.txt'
    shell:
        r"""
        tsinfer match-ancestors {input.samples} \
                -A {output} \
                --log-section tsinfer.inference \
                --num-threads {threads} \
                --progress \
                --verbosity \
                --ancestors {input.truncated_ancestors} 
        """
# --recombination-rate {params.re} \
               # --mismatch-ratio {params.mima} \
               

rule infer_trees:
    input:
        samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.samples',
        ancestors_trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.atrees',
    output:
        trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.trees',
    params:
        re = '1.1e-8',
        mima = '1'
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark:
            'benchmarks/{chromosome}_infer.benchmark.txt'
    shell:
        r"""
        # Check amount of memory (in kbytes) as seen by the job
        ulimit -v
        tsinfer match-samples {input.samples} \
                -A {input.ancestors_trees} \
                --num-threads {threads} \
                --progress \
                --log-section tsinfer.inference \
                --verbosity \
                -O {output}
        """
# --recombination-rate {params.re} \
#                 --mismatch-ratio {params.mima} \

rule simplify_trees:
    input: 
        trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.trees',
    output: 
        simplified_trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.simplified.trees',
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark:
            'benchmarks/{chromosome}_simplify.benchmark.txt'
    shell: 
        r"""
        ulimit -v
        python scripts/simplify.py {input} {output}
        """             
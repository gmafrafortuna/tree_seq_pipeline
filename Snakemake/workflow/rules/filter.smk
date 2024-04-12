if len(allFiles) != 1:
    rule compare_and_filter_vcfs:
        input: 
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.vcf.gz',
                    file=allFiles)]
        output:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz',
                    file=allFiles)]
        params:
            file_names = allFiles,
            path = f'{vcfdir}/{{chromosome}}/{{chromosome}}'
        run:
            import os
            import itertools
            import subprocess
            from glob import glob

            path_names = [f'{params.path}/{name}' for name in params.file_names]
            pairs = list(itertools.combinations(path_names, 2))
            
            for i, pair in enumerate(pairs):
                print(f'pair {i}')
                # for each pair determine the extension of the vcf files
                _0, ext0 = [os.path.splitext(file) for file in glob(f'{pair[0]}*') 
                                if os.path.splitext(file)[1] not in ['.csi', '.txt']][0]

                _1, ext1 = [os.path.splitext(file) for file in glob(f'{pair[1]}*') 
                                if os.path.splitext(file)[1] not in ['.csi', '.txt']][0]
                
                # Run vcftools to find if there are overlapping samples between pairs
                command = f"""
                        module load anaconda
                        conda activate bcftools
                        vcftools --gzvcf {os.path.join(pair[0] + ext0)} \
                                --gzdiff {os.path.join(pair[1] + ext1)} \
                                --diff-indv
                        """
                subprocess.run(command, shell=True, check=True, cwd=path)   

                # if vcftools find any overlapping samples, remove them from the second vcf file
                if os.path.exists(os.path.join(path, 'out.diff.indv_in_files')):
                    print(f'Removing samples from {pair[1]} and renaming {pair[0]}')
                    command = f"""
                            grep B out.diff.indv_in_files | cut -f1 > samples.remove
                            rm out.diff.indv_in_files
                            """
                    subprocess.run(command, shell=True, check=True, cwd=path)     
                    command = f"""
                              module load anaconda
                              conda activate bcftools
                              bcftools view -S ^samples.remove {os.path.join(pair[1] + '.vcf' + ext1)} \
                                            -Oz -o {os.path.join(pair[1] + new_ext)}
                              rm samples.remove
                            """
                    subprocess.run(command, shell=True, check=True, cwd=path)

                # if there are no overlapping samples, 
                # just rename the files that still have the original extension
                else:
                    print('No overlapping samples')
                    if 'filtered' not in _0:
                        print('Rename')
                        command = f"""
                                  mv {os.path.join(pair[0] + '.vcf' + ext0)} {os.path.join(pair[0] + new_ext)}
                                  """
                        subprocess.run(command, shell=True, check=True, cwd=path)           
                    if 'filtered' not in _1:
                        print('Rename')
                        command = f"""
                                  mv {os.path.join(pair[1] + '.vcf' + ext1)} {os.path.join(pair[1] + new_ext)}
                                  """
                        subprocess.run(command, shell=True, check=True, cwd=path)
                        
    # rule get_samples:
    #     input:
    #         f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz'
    #     output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
    #     conda: "bcftools"
    #     threads: 1
    #     resources: cpus=1, mem_mb=32000, time_min=30
    #     log: 'logs/Get_samples_{chromosome}_{file}.log'
    #     shell:
    #         """
    #         bcftools query -l {input} > {output}
    #         """

    # rule filter:
    #     input:
    #         [f'{vcfdir}' + x
    #             for x in expand('/{{chromosome}}/{{chromosome}}_{file}.txt',
    #                 file=allFiles)],
    #     output:
    #         temp([f'{vcfdir}' + x for x in
    #             expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.{ext}', 
    #             file=allFiles, ext=['gz', 'gz.csi'])])
    #     params:
    #         duplicated = '{chromosome}.ids',
    #         files_len  = '{chromosome}_len.txt',
    #         files_reo = '{chromosome}_reorderd.txt', 
    #     resources: cpus=1, mem_mb=64000, time_min=60
    #     run:
    #         import os

    #         for file in input:
    #             shell("wc -l {file} >> {params.files_len}")
            
    #         shell("sort -us -nr {params.files_len} | cut -d' ' -f 2 > {params.files_reo}")

    #         with open(params.files_reo, 'r') as file:
    #             files_reord = []
    #             for line in file:
    #                 line = line.strip()
    #                 files_reord.append(line)
            
    #         filtered = []

    #         for i in range(len(files_reord)):
    #             print(i)
    #             for j, file in enumerate(files_reord):
    #                 if file != files_reord[i] and j > i:
    #                     file1 = files_reord[i]
    #                     file2 = files_reord[j]
    #                     print(file1, file2)

    #                     shell("comm -12 <(sort {file1}) <(sort {file2}) > {params.duplicated}")

    #                     remove = params.duplicated
    #                     lremove = os.stat(f'{remove}').st_size

    #                     if lremove == 0:
    #                         print(f'Nothing to remove {lremove}')
    #                     else:
    #                         print(f'{lremove} samples to remove')
    #                         if file2 not in filtered:
    #                             filtered.append(file2)

    #                         file2 = file2.split('.')[0]
    #                         print(f'file {file2} filtered')
    #                         shell('bcftools view -S ^{params.duplicated} {file2}.vcf.gz -O z -o {file2}.filtered.vcf.gz')
    #                         #shell('gatk SelectVariants -R {params.refgen} -V {file2} -xl-sn -O {file2}.filtered.vcf.gz')
    #                         shell('bcftools index -f {file2}.filtered.vcf.gz')
    #                         shell('rm {params.duplicated} {file2}.vcf.gz {file2}.vcf.gz.csi')

    #         for file in files_reord:
    #             if file not in filtered:
    #                 file1 = file.split('.')[0] # file does not have extensions
    #                 print(f'file {file1} renamed')   
    #                 shell('mv {file1}.vcf.gz {file1}.filtered.vcf.gz')
    #                 shell('mv {file1}.vcf.gz.csi {file1}.filtered.vcf.gz.csi')
    #             #shell('scripts/create_symlinks.sh {file1}')
    #         shell('rm {params}')

else:
    rule rename:
        input:
            vcf = [[f'{vcfdir}' + x
                    for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                        suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
                    [f'{vcfdir}' + x
                        for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                            suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []],
            idx = [[f'{vcfdir}' + x
                    for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz.csi',
                        suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
                    [f'{vcfdir}' + x
                        for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz.csi',
                            suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []],
        output:
            vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.filtered.vcf.gz',
            idx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.filtered.vcf.gz.csi'
        conda: "bcftools"
        log: 'logs/Rename_{chromosome}.log'
        resources: cpus=1, mem_mb=32000, time_min=30
        shell:
            """
            if [ -h {input} ]; then
                ln -s $( realpath {input.vcf} ) {output.vcf}
                ln -s $( realpath {input.idx} ) {output.idx}
            else
                ln -s {input.vcf} {output.vcf}
                ln -s {input.idx} {output.idx}
            fi
            """


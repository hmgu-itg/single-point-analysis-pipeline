"""
Snakefile 5. Run COJO

Example
-------
# To run cojo on all detected peaks
$ snakemake --cores 10 --keep-going --use-singularity --snakefile workflow/rules/cojo.smk run_all_cojo

"""
include: "1. single-cohort.smk"

import pandas as pd

def run_all_cojo_input(w):
    peaklist = checkpoints.detect_peaks.get(cohort=w.cohort, group=w.group, phenotype=w.phenotype).output[0]
    try:
        peaklist = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
    except pd.errors.EmptyDataError:
        return []
    peaks = [f'output/{w.cohort}/{group}/{phenotype}/cojo/{chrom}.{start}.{end}.jma.cojo' for _, (group, phenotype, chrom, start, end) in peaklist.iterrows()]
    return peaks

rule run_all_cojo:
    input:
        run_all_cojo_input
    output:
        f"{OUTPUT_PATH}/all.cojo.jma.csv.gz"
    run:
        if len(input)==0: # When no signif signal
            empty = pd.DataFrame(columns = ["group", "phenotype", "peak", "Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r", "cojo_tophit"])
            empty.to_csv(output[0], index = False, header = True, compression = 'gzip')
            return

        concat_list = list()
        for f in input:
            string = re.sub(r'.*\/cojo\/', '', f)
            string = re.sub(r'\/.*', '', string)
            peak = f.split('/')[-1].replace(f'{string}.', '').replace('.jma.cojo', '')
            df = pd.read_csv(f, sep = '\t', header = 0)
            min_p = df['p'].min()
            df.loc[df['p'] == min_p, 'cojo_tophit'] = True
            df.insert(0, 'group', wildcards.group)
            df.insert(1, 'phenotype', wildcards.phenotype)
            df.insert(2, 'peak', peak)
            concat_list.append(df)
        pd.concat(concat_list).to_csv(output[0], index = False, header = True, compression = 'gzip')


# rule get_samplesize:
#     '''
#     This file is used when creating the cojofile.    
#     Warning:
#     The awk command was only tested in mawk (1.3.4 20200120). 
#     This awk command may not work as intended using other type of awk (original awk, gawk)!
#     '''
#     input:
#         BFILE_INPUTS
#     params:
#         input=BFILE,
#         output='output/{cohort}/{cohort}.miss'
#     threads: workflow.cores
#     output:
#         lmiss=temp('output/{cohort}/{cohort}.miss.lmiss'),
#         imiss=temp('output/{cohort}/{cohort}.miss.imiss'),
#         samplesize='output/{cohort}/{cohort}.samplesize'
#     shell: """
#         plink \
#           --bfile {params.input} \
#           --memory {resources.mem_mb} \
#           --threads {threads} \
#           --missing \
#           --out {params.output} \
#           --silent

#         awk -F '[[:space:]]+' '{{if(NR!=1){{print $3, $5-$4}}}}' {output.lmiss} > {output.samplesize}
#         """


rule get_samplesize:
    '''
    This file is used when creating the cojofile.    
    Warning:
    The awk command was only tested in mawk (1.3.4 20200120). 
    This awk command may not work as intended using other type of awk (original awk, gawk)!
    '''
    input:
        lmiss=LMISS
    output:
        samplesize=f"{OUTPUT_PATH}/samplesize.txt"
    shell: """
        awk -F '[[:space:]]+' '{{if(NR!=1){{print $3, $5-$4}}}}' {input.lmiss} > {output.samplesize}
    """

rule make_cojofile:
    '''
    Creates cojofile from bfile and assoc file.
    bfile is only used to calculate the "sample size" value for each variant (The "N" column of cojofile).
    The input assoc file makes up bulk of the cojofile. 
    
    Reference
    --------
    1. https://yanglab.westlake.edu.cn/software/gcta/#COJO
    2. https://yanglab.westlake.edu.cn/software/gcta/#MLMA
    '''
    input:
        assoc=rules.gcta.output.mlma_bgz,
        samplesize=rules.get_samplesize.output.samplesize
    # params:
    #     bfile=BFILE,
    #     prefix=f"{OUTPUT_PATH}/cojofile/cojofile"
    output:
        f"{OUTPUT_PATH}/cojofile/cojofile.ma.gz"
    run:
        samplesize = pd.read_csv(input.samplesize,
                                sep = ' ',
                                header = None,
                                names = ['SNP', 'N'])
        assoc = pd.read_csv(input.assoc,
                            sep = '\t',
                            usecols = ['SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p'])
        cojofile = assoc.merge(samplesize)
        assert assoc.shape[0] == cojofile.shape[0], 'association file and cojofile row count do not match!'
        cojofile.rename(columns = {'Freq': 'freq'}, inplace=True)

        cojofile = cojofile[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']] # Just to be absolutely sure
        cojofile.to_csv(output[0],
                        sep = ' ',
                        index = False,
                        compression = 'gzip')



rule cojo:
    input:
        bfiles=BFILE_INPUTS,
        cojofile=rules.make_cojofile.output,
        assoc=rules.gcta.output.mlma_bgz
    params:
        bfile=BFILE,
        threshold=config['p-value'],
        prefix=f"{OUTPUT_PATH}/cojo/{{chrom}}.{{start}}.{{end}}"
    output:
        multiext(f"{OUTPUT_PATH}/cojo/{{chrom}}.{{start}}.{{end}}", 
            ".jma.cojo",
            ".cma.cojo",
            ".ldr.cojo")
    log:
        f"{OUTPUT_PATH}/cojo/{{chrom}}.{{start}}.{{end}}.cojo.log"
    shell:
        """
        workflow/scripts/run_cojo.sh \
            {params.bfile} \
            {input.cojofile} \
            {input.assoc} \
            {params.threshold} \
            {wildcards.chrom} \
            {wildcards.start} \
            {wildcards.end} \
            {params.prefix} \
            {resources.mem_mb} 2>&1 > {log}
        """

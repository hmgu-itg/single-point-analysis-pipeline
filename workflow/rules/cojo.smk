"""
Snakefile 5. Run COJO

Example
-------
# To run cojo on all detected peaks
$ snakemake --cores 10 --keep-going --use-singularity --snakefile workflow/rules/cojo.smk

"""
import re

import pandas as pd

include: "read-config.smk"
include: "meta-analysis.smk"

container: config['container']

peaklist = expand("output/meta-analysis/peaks/peaklist/all.{group}.peaklist", group=config['group'])[0]

all_peaks = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
peaks = [f'{row.group}.{row.phenotype}/{row.group}.{row.phenotype}.{row.chrom}.{row.start}.{row.end}' for _, row in all_peaks.iterrows()]

rule all_cojo:
    input:
        expand("output/meta-analysis/cojo/{peak}.jma.cojo", peak=peaks)
    output:
        "output/meta-analysis/cojo/all.jma.cojo.csv.gz"
    run:
        concat_list = list()
        for f in input:
            string = re.sub(r'.*\/cojo\/', '', f)
            string = re.sub(r'\/.*', '', string)
            group, phenotype = string.split('.')
            peak = f.split('/')[-1].replace(f'{string}.', '').replace('.jma.cojo', '')
            df = pd.read_csv(f, sep = '\t', header = 0)
            min_p = df['p'].min()
            df['is cojo tophit'] = False
            df.loc[df['p'] == min_p, 'is cojo tophit'] = True
            df.insert(0, 'peak', peak)
            df.insert(0, 'phenotype', phenotype)
            df.insert(0, 'group', group)
            concat_list.append(df)
        pd.concat(concat_list).to_csv(output[0], index = False, header = True, compression = 'gzip')


rule make_cojofile:
    '''
    Creates cojofile from bfile and metal.
    bfile is only used to calculate the "sample size" value for each variant (The "N" column of cojofile).
    The input metal file makes up bulk of the cojofile. 
    

    Reference
    --------
    1. https://yanglab.westlake.edu.cn/software/gcta/#COJO
    2. https://genome.sph.umich.edu/wiki/METAL_Documentation
    '''
    input:
        metal=rules.filter_metal.output.gz,
        bfiles=rules.phenotype_mac_filter.output.bfile
    params:
        group="{group}",
        phenotype="{phenotype}",
        bfile=rules.phenotype_mac_filter.params.out,
        prefix="output/meta-analysis/cojofile/{group}.{phenotype}"
    output:
        "output/meta-analysis/cojofile/{group}.{phenotype}.ma.gz"
    log:
        "output/meta-analysis/cojofile/{group}.{phenotype}.log"
    shell:
        """
        workflow/scripts/make_cojofile.sh \
            {params.group} \
            {params.phenotype} \
            {params.bfile} \
            {input.metal} \
            {params.prefix} 2>&1 > {log}
        """


rule cojo:
    input:
        bfiles=rules.phenotype_mac_filter.output.bfile,
        cojofile=rules.make_cojofile.output,
        metal=rules.filter_metal.output.gz
    params:
        group="{group}",
        phenotype="{phenotype}",
        bfile="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        threshold=config['QC_thresholds']['p-value'],
        chrom="{chrom}",
        start="{start}",
        end="{end}",
        prefix="output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}"
    output:
        multiext("output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}", 
            ".jma.cojo",
            ".cma.cojo",
            ".ldr.cojo")
    log:
        "output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}.cojo.log"
    shell:
        """
        workflow/scripts/run_cojo.sh \
            {params.bfile} \
            {input.cojofile} \
            {input.metal} \
            {params.threshold} \
            {params.group} \
            {params.phenotype} \
            {params.chrom} \
            {params.start} \
            {params.end} \
            {params.prefix} 2>&1 > {log}
        """

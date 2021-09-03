import pandas as pd
include: "read-config.smk"
container: config['container']

peaklist = expand("output/meta-analysis/peaks/peaklist/all.{group}.peaklist", group=config['group'])[0]

all_peaks = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
peaks = [f'{row.group}.{row.phenotype}/{row.group}.{row.phenotype}.{row.chrom}.{row.start}.{row.end}' for _, row in all_peaks.iterrows()]

rule all_cojo:
    input:
        expand("output/meta-analysis/cojo/{peak}.jma.cojo", peak=peaks)



rule make_cojofile:
    input:
        metal="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz",
        metal_index="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz.tbi",
        bfiles=multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam')
    params:
        group="{group}",
        phenotype="{phenotype}",
        bfile="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        prefix="output/meta-analysis/cojofile/{group}.{phenotype}"
    output:
        "output/meta-analysis/cojofile/{group}.{phenotype}.miss500.ma.gz"
    log:
        "output/meta-analysis/cojofile/{group}.{phenotype}.miss500.ma.log"
    shell:
        """
        workflow/scripts/make_cojofile.sh \
            {params.group} \
            {params.phenotype} \
            {params.bfile} \
            {input.metal} \
            {params.prefix} 2>&1 > {log}
        """


rule run_cojo:
    input:
        bfiles=multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam'),
        cojofile="output/meta-analysis/cojofile/{group}.{phenotype}.miss500.ma.gz",
        metal="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz",
    params:
        group="{group}",
        phenotype="{phenotype}",
        bfile="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        chrom="{chrom}",
        start="{start}",
        end="{end}",
        prefix="output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}"
    output:
        multiext("output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}", 
            ".jma.cojo",
            ".cma.cojo",
            ".ldr.cojo",
            ".qassoc")
    log:
        "output/meta-analysis/cojo/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}.snakemake.log"
    shell:
        """
        workflow/scripts/run_cojo.sh \
            {params.bfile} \
            {input.cojofile} \
            {input.metal} \
            {params.group} \
            {params.phenotype} \
            {params.chrom} \
            {params.start} \
            {params.end} \
            {params.prefix} 2>&1 > {log}
        """
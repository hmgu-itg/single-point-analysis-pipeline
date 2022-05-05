"""
Snakefile 3. Detect peaks


Example
-------
$ snakemake --cores 10 --keep-going --use-singularity --snakefile workflow/rules/detect-peaks.smk  detected_all_peaks
$ snakemake --profile slurm --keep-going --use-singularity --snakefile workflow/rules/detect-peaks.smk detected_all_peaks
"""
include: "read-config.smk"
include: 'meta-analysis.smk'

rule detect_all_peaks:
    input:
        expand("output/meta-analysis/peaks/peaklist/all.{group}.peaklist", group = config['group'])

rule collate_all_peaks:
    input:
        expand("output/meta-analysis/peaks/peaklist/{group}.{phenotype}.peaklist", 
                phenotype = config['phenotypes'], allow_missing=True)
    output:
        "output/meta-analysis/peaks/peaklist/all.{group}.peaklist"
    resources:
        mem_mb=500,
        time="0:30:0",
        cpus_per_task=1
    run:
        all_df = list()
        for f in input:
            try:
                df = pd.read_csv(f, sep = '\t', header = None)
            except pd.errors.EmptyDataError:
                continue
            all_df.append(df)
        all_df = pd.concat(all_df)
        all_df.sort_values([0, 1, 2, 3, 4], inplace = True)
        all_df.to_csv(output[0], sep = '\t', header = False, index = False)


rule detect_peaks:
    input:
        rules.filter_metal.output.gz
        # "output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz"
    params:
        span=config['peakplotter']['span'],
        signif=config['QC_thresholds']['p-value'],
        group="{group}",
        phenotype="{phenotype}"
    output:
        "output/meta-analysis/peaks/peaklist/{group}.{phenotype}.peaklist"
    log:
        "output/meta-analysis/peaks/peaklist/{group}.{phenotype}.log"
    singularity: "library://hmgu-itg/default/peakplotter:0.4.3"
    resources:
        mem_mb=500,
        time="0:30:0",
        cpus_per_task=1
    shell:
        "python3 workflow/scripts/collect_peaks.py {input} {params.span} {params.signif} {params.group} {params.phenotype} {output} 2>&1 > {log}"

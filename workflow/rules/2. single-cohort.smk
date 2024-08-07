"""
Snakefile 1.

Run association and create QQ-plot and PeakPlotter plot.
"""
include: "0. read-config.smk"

        
rule gcta:
    input:
        phenotype=PHENOTYPE_FILE,
        bfile=BFILE,
        grm=multiext(GRM, '.grm.bin', '.grm.id', '.grm.N.bin')
    params:
        out="output/{cohort}/{group}/{phenotype}/gcta/gcta",
        bfile="output/{cohort}/bfile/{cohort}",
        grm=GRM
    output:
        pheno="output/{cohort}/{group}/{phenotype}/gcta/gcta.pheno",
        mlma=temp("output/{cohort}/{group}/{phenotype}/gcta/gcta.mlma"),
        mlma_bgz="output/{cohort}/{group}/{phenotype}/gcta/gcta.mlma.gz",
        mlma_bgz_tbi="output/{cohort}/{group}/{phenotype}/gcta/gcta.mlma.gz.tbi"
    threads: workflow.cores
    log: "output/{cohort}/{group}/{phenotype}/gcta/gcta.mlma.log"
    shell:
        """
        awk '{{OFS=\"\\t\"}}{{print $1,$1,$2}}' {input.phenotype} > {output.pheno} 2> {log}
        gcta64 --mlma \
                --bfile {params.bfile} \
                --grm {params.grm} \
                --pheno {output.pheno} \
                --out {params.out} \
                --threads {threads} 2>&1 >> {log}
        bgzip -c {output.mlma} > {output.mlma_bgz} 2>> {log}
        tabix --skip-lines 1 --sequence 1 --begin 3 --end 3 {output.mlma_bgz} 2>&1 >> {log}
        
        rm {params.out}.log 2>&1 >> {log}
        """


rule manqq_gcta:
    input: rules.gcta.output.mlma_bgz
    params:
        prefix="output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}",
        filter="{filter}"
    resources:
        rate_limit=1
    output: 
        "output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}.run_conf",
        "output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}.qq.png",
        "output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}.lambda.txt"
    log:
        out="output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}.o",
        err="output/{cohort}/{group}/{phenotype}/manqq/manqq.{filter}.e"
    singularity: config['container']['manqq']
    shell:
        """
        manqq_cli \
          --chr-col Chr \
          --pval-col p \
          --pos-col bp \
          --a1 A1 \
          --a2 A2 \
          --build 38 \
          --image png \
          --af-col Freq \
          --qq-title {wildcards.phenotype} \
          --manh-title {wildcards.phenotype} \
          --maf-filter {params.filter} \
          {input} \
          {params.prefix} > {log.out} 2> {log.err}
        """


checkpoint detect_peaks:
    input:
        rules.gcta.output.mlma_bgz
    params:
        span=config['peakplotter']['span'],
        signif=config['p-value']
    output:
        "output/{cohort}/{group}/{phenotype}/peaklist"
    log:
        "output/{cohort}/{group}/{phenotype}/peaklist.log"
    singularity: config['container']['peakplotter']
    shell:
        "python3 workflow/scripts/collect_peaks.py {input} {params.span} {params.signif} {wildcards.group} {wildcards.phenotype} {output} 2>&1 > {log}"


rule plotpeak:
    input:
        assoc=rules.gcta.output.mlma_bgz,
        bfile=BFILE_INPUTS
    params:
        bfile=BFILE,
        outdir="output/{cohort}/{group}/{phenotype}/peaks",
        chrom="{chrom}",
        start="{start}",
        end="{end}",
        vep_ld=config['peakplotter']['vep_ld']
    singularity: config['container']['peakplotter']
    resources:
        rate_limit=1
    output:
        multiext("output/{cohort}/{group}/{phenotype}/peaks/{chrom}.{start}.{end}.500kb", '.html', '.csv')
    shell:
        """
        peakplotter-region  \
          --chr-col Chr \
          --pos-col bp \
          --rs-col SNP \
          --pval-col p \
          --a1-col A1 \
          --a2-col A2 \
          --maf-col Freq \
          --build 38 \
          --assoc-file {input.assoc} \
          --bfiles {params.bfile} \
          --out {params.outdir} \
          --chrom {wildcards.chrom} \
          --start {wildcards.start} \
          --end {wildcards.end} \
          --vep-ld {params.vep_ld} \
          --debug
        """


def plot_all_peaks_input(w):
    peaklist = checkpoints.detect_peaks.get(cohort=w.cohort, group=w.group, phenotype=w.phenotype).output[0]
    try:
        peaklist = pd.read_csv(peaklist, sep = '\t', header = None)
    except pd.errors.EmptyDataError:
        return []
    return [rules.plotpeak.output[0].format(cohort=w.cohort, group=group, phenotype=phenotype, chrom=chrom, start=start, end=end)
            for _, (group, phenotype, chrom, start, end) in peaklist.iterrows()]

rule plot_all_peaks:
    input:
        plot_all_peaks_input
    output:
        touch("output/{cohort}/{group}/{phenotype}/peaks/.done")
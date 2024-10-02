"""
Snakefile 1.

Run association and create QQ-plot and PeakPlotter plot.
"""
configfile: "config.yaml"
container: config['container']['all']


rule gcta:
    input:
        phenotype=PHENOTYPE_FILE,
        bfile=multiext(config['bfile'], '.bed', '.bim', '.fam'),
        grm=multiext(config['grm'], '.grm.bin', '.grm.id', '.grm.N.bin')
    params:
        out="{output}/gcta/gcta",
        bfile=config['bfile'],
        grm=config['grm']
    output:
        pheno="{output}/gcta/gcta.pheno",
        mlma=temp("{output}/gcta/gcta.mlma"),
        mlma_bgz="{output}/gcta/gcta.mlma.gz",
        mlma_bgz_tbi="{output}/gcta/gcta.mlma.gz.tbi"
    threads: workflow.cores
    log: "{output}/gcta/gcta.mlma.log"
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
        prefix="{output}/manqq/manqq.{filter}"
    resources:
        rate_limit=1
    output: 
        "{output}/manqq/manqq.{filter}.run_conf",
        "{output}/manqq/manqq.{filter}.qq.png",
        "{output}/manqq/manqq.{filter}.lambda.txt"
    log:
        out="{output}/manqq/manqq.{filter}.o",
        err="{output}/manqq/manqq.{filter}.e"
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
          --maf-filter {wildcards.filter} \
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
        "{output}/peaklist"
    log:
        "{output}/peaklist.log"
    singularity: config['container']['peakplotter']
    shell:
        "python3 workflow/scripts/collect_peaks.py {input} {params.span} {params.signif} {wildcards.group} {wildcards.phenotype} {output} 2>&1 > {log}"


rule plotpeak:
    input:
        assoc=rules.gcta.output.mlma_bgz,
        bfile=multiext(config['bfile'], '.bed', '.bim', '.fam')
    params:
        bfile=config['bfile'],
        outdir="{output}/peaks",
        chrom="{chrom}",
        start="{start}",
        end="{end}",
        vep_ld=config['peakplotter']['vep_ld']
    singularity: config['container']['peakplotter']
    resources:
        rate_limit=1
    output:
        multiext("{output}/peaks/{chrom}.{start}.{end}.500kb", '.html', '.csv')
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
    peaklist = checkpoints.detect_peaks.get().output[0]
    try:
        peaklist = pd.read_csv(peaklist, sep = '\t', header = None)
    except pd.errors.EmptyDataError:
        return []
    return [rules.plotpeak.output[0].format(chrom=chrom, start=start, end=end)
            for _, (group, phenotype, chrom, start, end) in peaklist.iterrows()]

rule plot_all_peaks:
    input:
        plot_all_peaks_input
    output:
        touch("{output}/peaks/.done")

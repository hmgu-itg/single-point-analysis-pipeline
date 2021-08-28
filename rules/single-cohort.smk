configfile: "config.yaml"


rule all:
    expand("output/gcta/{cohort}/{cohort}.{group}.{phenotype}", cohort=config['cohorts'], group=config['group'], phenotype=config['phenotypes'])

rule gcta:
    input:
        phenotype=config['input'],
        bfile=multiext("output/bfile/{cohort}", '.bed', '.bim', '.fam'),
        grm=lambda w: multiext(config["cohorts"][w.cohort]['grm'], '.grm.bin', '.grm.id', '.grm.N.bin')
    params:
        outprefix="output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}",
        bfile="output/bfile/{cohort}",
        grm=lambda w: config["cohorts"][w.cohort]['grm']
    output:
        mlma_bgz="output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz",
        mlma_bgz_tbi="output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz.tbi"
    threads: 10
    log: "output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.log"
    shell:
        """
        scripts/run_gcta.sh \
            {input.phenotype} \
            {params.outprefix} \
            {params.bfile} \
            {params.grm} \
            {threads} 2>&1 > {log}
        rm {params.outprefix}.log
        """


rule peak_gcta:
    input:
        mlma="output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz",
        bfile=multiext("output/bfile/{cohort}", '.bed', '.bim', '.fam')
    singularity: "library://hmgu-itg/default/peakplotter"
    params:
        bfile="output/bfile/{cohort}"
    resources:
        rate_limit=1
    output: 
        directory=directory("output/single-cohort/peaks/{cohort}/{cohort}.{panel}.{protein}"),
        done="output/single-cohort/peaks/{cohort}/{cohort}.{panel}.{protein}/done"
    log: "output/single-cohort/peaks/{cohort}/{cohort}.{panel}.{protein}"
    shell:
        """
        peakplotter  \
          --bfiles {params.bfile} \
          --chr-col Chr \
          --pos-col bp \
          --rs-col SNP \
          --pval-col p \
          --a1-col A1 \
          --a2-col A2 \
          --maf-col Freq \
          --build 38 \
          --overwrite \
          --assoc-file {input.mlma} \
          --out {output.directory} 2>&1 > {log}
        touch {output.done}
        """
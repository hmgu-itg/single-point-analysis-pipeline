"""
Snakefile 1.

Run GCTA-MLMA and ManQQ

Example
-------
# To run GCTA and ManQQ on all
$ snakemake --cores 10 --restart-times 3 --keep-going --use-singularity --snakefile workflow/rules/single-cohort.smk 

# To run a single GCTA run
$ snakemake \
    --cores 1 \
    --restart-times 3 \
    --keep-going \
    --use-singularity \
    --snakefile workflow/rules/single-cohort.smk \
    output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz

# To run a single ManQQ run
$ snakemake \
    --cores 1 \
    --restart-times 3 \
    --keep-going \
    --use-singularity \
    --snakefile workflow/rules/single-cohort.smk \
    output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.qq.png

"""
include: "read-config.smk"

rule all2:
    input:
        expand("output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz", cohort=config['cohorts'], group=config['group'], phenotype=config['phenotypes']),
        expand("output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.qq.png", 
               cohort=config['cohorts'], group=config['group'], phenotype=config['phenotypes'], filter=[0.0, 0.001])
        
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
    threads: 20
    resources:
        cpus_per_task=20,
        mem_per_cpu="10G"
    log: "output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.log"
    shell:
        """
        pheno_file={params.outprefix}.pheno
        mlma={params.outprefix}.mlma
        mlma_bgz={params.outprefix}.mlma.gz

        awk '{{OFS=\"\\t\"}}NR!=1{{print $1,$1,$3}}' {input.phenotype} > $pheno_file 2> {log}

        gcta64 --mlma \
                --bfile {params.bfile} \
                --grm {params.grm} \
                --pheno $pheno_file \
                --out {params.outprefix} \
                --threads {threads} 2>&1 >> {log}
        bgzip -c $mlma > $mlma_bgz 2>> {log}
        tabix --skip-lines 1 --sequence 1 --begin 3 --end 3 $mlma_bgz 2>&1 >> {log}

        rm $pheno_file $mlma {params.outprefix}.log 2>&1 >> {log}
        """


rule manqq_gcta:
    input: "output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz"
    params:
        prefix="output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}",
        filter="{filter}"
    resources:
        cpus_per_task=1,
        mem_per_cpu="3G",
        rate_limit=1
    output: 
        "output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.run_conf",
        "output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.qq.png",
        "output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.lambda.txt"
    log:
        out="output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.o",
        err="output/single-cohort/manqq/{cohort}/{cohort}.{group}.{phenotype}.manqq.{filter}.e"
    container: "library://hmgu-itg/default/manqq:0.2.3"
    shell:
        """
        run_manqq.R \
          --chr-col Chr \
          --pval-col p \
          --pos-col bp \
          --a1 A1 \
          --a2 A2 \
          --build 38 \
          --image png \
          --af-col Freq \
          --no-man \
          --maf-filter {params.filter} \
          {input} \
          {params.prefix} > {log.out} 2> {log.err}
        """

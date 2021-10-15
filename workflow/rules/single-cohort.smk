"""
Snakefile 1.


"""
include: "read-config.smk"

rule all2:
    input:
        expand("output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz", cohort=config['cohorts'], group=config['group'], phenotype=config['phenotypes'])

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

# configfile: "config.yaml"
include: "read-config.smk"
# conda: "conda/environment.yaml"
# snakemake --profile slurm --use-conda --keep-going --show-failed-log --quiet --snakefile rules/single-cohort.smk --conda-frontend conda

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
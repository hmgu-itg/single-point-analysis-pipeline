# configfile: "config.yaml"
include: "read-config.smk"

rule all:
    input:
        expand("output/bfile/{cohort}.bed", cohort=config['cohorts'])
        'output/bfile/all.mac.txt'


rule hwe:
    input: lambda w: multiext(config['cohorts'][w.cohort]['bfile'], '.bed', '.bim', '.fam')
    params:
        input=lambda w: config['cohorts'][w.cohort]['bfile'],
        output="output/bfile/hwe/{cohort}"
    output: multiext("output/bfile/hwe/{cohort}", '.hwe', '.log', '.nosex')
    shell:
        """
        plink \
          --bfile {params.input} \
          --hardy \
          --out {params.output}
        """

rule exclude_hwe:
    input: "output/bfile/hwe/{cohort}.hwe"
    params: config['QC_thresholds']['HWE']
    output: "output/bfile/hwe/{cohort}.exclude.txt"
    shell:
        """
        tail -n+2 {input} | tr -s ' ' | awk '{{if ($NF<={params}) {{print $2}}}}' > {output}
        """

rule missingness:
    input: lambda w: multiext(config['cohorts'][w.cohort]['bfile'], '.bed', '.bim', '.fam')
    params:
        input=lambda w: config['cohorts'][w.cohort]['bfile'],
        output="output/bfile/missingness/{cohort}"
    output: multiext("output/bfile/missingness/{cohort}", '.lmiss', '.imiss')
    shell:
        """
        plink \
          --bfile {params.input} \
          --missing \
          --out {params.output}
        """

rule exclude_missingness:
    input: "output/bfile/missingness/{cohort}.lmiss"
    params: config['QC_thresholds']['missingness']
    output: "output/bfile/missingness/{cohort}.exclude.txt"
    shell:
        """
        tail -n+2 {input} | tr -s ' ' | awk '{{if ($NF>={params}) {{print $2}}}}' > {output}
        """


rule collect_exclude_list:
    input:
        "output/bfile/missingness/{cohort}.exclude.txt",
        "output/bfile/hwe/{cohort}.exclude.txt",
    output: "output/bfile/{cohort}.exclude.all.txt"
    shell:
        "cat {input} > {output}"


rule filter_bfile:
    input: 
        bfile=lambda w: multiext(config['cohorts'][w.cohort]['bfile'], '.bed', '.bim', '.fam'),
        exclude_list="output/bfile/{cohort}.exclude.all.txt"
    params:
        bfile=lambda w: config['cohorts'][w.cohort]['bfile'],
        out="output/bfile/{cohort}"
    output: multiext("output/bfile/{cohort}", '.bed', '.bim', '.fam')
    threads: 10
    shell:
        """
        plink \
          --bfile {params.bfile} \
          --exclude {input.exclude_list} \
          --make-bed \
          --out {params.out} \
          --threads {threads}
        awk -F$'\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1,\"chr\"$1\":\"$4,$3,$4,$5,$6}}' {params.out}.bim | sponge {params.out}.bim
        """


rule mac:
    input: config['input']
    params:
        id="{cohort}.{group}.{phenotype}",
        threshold=config['QC_thresholds']['MAC']
    output: temp('output/bfile/mac/{cohort}.{group}.{phenotype}.txt')
    run:
        data = pd.read_csv(input[0], sep = '\t')
        n = data.shape[0]
        max_mac = n * 2
        threshold = params.threshold / max_mac
        pd.DataFrame([[params.id, n, max_mac, threshold]]).to_csv(output[0], index = False, header = False)


rule mac_group:
    input: lambda w: expand('output/bfile/mac/{cohort}.{group}.{phenotype}.txt', cohort=config['cohorts'], group=w.group, phenotype=config['phenotypes'][w.group])
    output: temp("output/bfile/mac/{group}.txt")
    shell: "cat {input} > {output}"
    

rule mac_all:
    input: expand('output/bfile/mac/{group}.txt', group=config['phenotypes'])
    output: 'output/bfile/all.mac.txt'
    shell:
        "cat {input} > {output}"





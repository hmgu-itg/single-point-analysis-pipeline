import pandas as pd

include: "0. read-config.smk"
container: config['container']['all']

rule step1:
    input:
        multiext(f"output/{COHORT}/bfile/{COHORT}", '.bed', '.bim', '.fam')


rule hwe:
    input: multiext(BFILE, '.bed', '.bim', '.fam')
    params:
        input=BFILE,
        output="output/{cohort}/bfile/{cohort}"
    threads: workflow.cores * 0.5
    output: multiext("output/{cohort}/bfile/{cohort}", '.hwe', '.log', '.nosex')
    shell:
        """
        plink \
          --bfile {params.input} \
          --memory {resources.mem_mb} \
          --threads {threads} \
          --hardy \
          --out {params.output}
        """

rule exclude_hwe:
    input: "output/{cohort}/bfile/{cohort}.hwe"
    params: config['QC_thresholds']['HWE']
    output: "output/{cohort}/bfile/{cohort}.hwe.exclude.txt"
    shell:
        """
        tail -n+2 {input} | tr -s ' ' | awk '{{if ($NF<={params}) {{print $2}}}}' > {output}
        """

rule missingness:
    input: multiext(BFILE, '.bed', '.bim', '.fam')
    params:
        input=BFILE,
        output="output/{cohort}/bfile/{cohort}"
    threads: workflow.cores * 0.5
    output:
        lmiss="output/{cohort}/bfile/{cohort}.lmiss",
        imiss="output/{cohort}/bfile/{cohort}.imiss"
    shell:
        """
        plink \
          --bfile {params.input} \
          --memory {resources.mem_mb} \
          --threads {threads} \
          --missing \
          --out {params.output}
        """

rule exclude_missingness:
    input: rules.missingness.output.lmiss
    params: config['QC_thresholds']['missingness']
    output: "output/{cohort}/bfile/{cohort}.missingness.exclude.txt"
    shell:
        """
        tail -n+2 {input} | tr -s ' ' | awk '{{if ($NF>={params}) {{print $2}}}}' > {output}
        """


rule collect_exclude_list:
    input:
        "output/{cohort}/bfile/{cohort}.hwe.exclude.txt",
        "output/{cohort}/bfile/{cohort}.missingness.exclude.txt",
    output: "output/{cohort}/bfile/{cohort}.all.exclude-list.txt"
    shell:
        "cat {input} > {output}"


rule filter_bfile:
    """
    Single-cohort bfiles are filtered on HWE and Missingness.
    """
    input: 
        bfile=multiext(BFILE, '.bed', '.bim', '.fam'),
        exclude_list="output/{cohort}/bfile/{cohort}.all.exclude-list.txt"
    params:
        bfile=BFILE,
        out="output/{cohort}/bfile/{cohort}"
    threads: workflow.cores
    output: multiext("output/{cohort}/bfile/{cohort}", '.bed', '.bim', '.fam')
    shell:
        """
        plink \
          --bfile {params.bfile} \
          --exclude {input.exclude_list} \
          --make-bed \
          --out {params.out} \
          --memory {resources.mem_mb} \
          --threads {threads}
        awk -F$'\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1,\"chr\"$1\":\"$4,$3,$4,$5,$6}}' {params.out}.bim | sponge {params.out}.bim
        """


rule get_samplesize:
    '''
    This file is used later when creating the cojofile.
    
    Warning:
    The awk command was only tested in mawk (1.3.4 20200120). 
    This awk command may not work as intended using other type of awk (original awk, gawk)!
    '''
    input:
        rules.filter_bfile.output
    params:
        input='output/{cohort}/bfile/{cohort}',
        output='output/{cohort}/bfile/{cohort}.miss'
    threads: workflow.cores
    output:
        lmiss=temp('output/{cohort}/bfile/{cohort}.miss.lmiss'),
        imiss=temp('output/{cohort}/bfile/{cohort}.miss.imiss'),
        samplesize='output/{cohort}/bfile/{cohort}.samplesize'
    shell: """
        plink \
          --bfile {params.input} \
          --memory {resources.mem_mb} \
          --threads {threads} \
          --missing \
          --out {params.output} \
          --silent

        awk -F '[[:space:]]+' '{{if(NR!=1){{print $3, $5-$4}}}}' {output.lmiss} > {output.samplesize}
        """




# rule mac:
#     input: PHENOTYPE_FILE
#     params:
#         phenotype=PHENOTYPE,
#         threshold=config['QC_thresholds']['MAC']
#     output: 'output/{cohort}/{group}/{phenotype}/{phenotype}.mac.txt'
#     run:
#         data = pd.read_csv(input[0], sep = '\t')
#         n = data.shape[0]
#         max_mac = n * 2
#         threshold = params.threshold / max_mac
#         pd.DataFrame([[params.phenotype, n, max_mac, threshold]]).to_csv(output[0], sep = ' ', index = False, header = False)

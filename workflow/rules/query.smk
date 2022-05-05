"""
Snakefile 6. Query GWAS Catalog and Ensembl VEP


"""
import json
import requests
import pandas as pd
include: "read-config.smk"


# chunk_list = [f.replace('.txt', '.csv') for f in pd.read_csv("output/meta-analysis/query/vep_chunks.txt", header = None)[0]]
peaklist = f"output/meta-analysis/peaks/peaklist/all.{config['group']}.peaklist"
peaklist = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
gwas_query_list = [f'output/meta-analysis/query/gwas/{r.group}.{r.phenotype}.{r.chrom}.{r.start}.{r.end}.csv' for _, r in peaklist.iterrows()]


rule query_all_gwas:
    input:
        gwas_query_list
    output:
        "output/meta-analysis/query/gwas.filter.yaml"
    run:
        concat_list = [pd.read_csv(f) for f in input]
        out = pd.concat(concat_list)

        with open(output[0], 'w') as f:
            for pheno in out['phenotype'].drop_duplicates():
                f.write(f"{pheno}: ['{pheno}']\n")

rule query_gwas:
    output:
        "output/meta-analysis/query/gwas/{group}.{phenotype}.{chrom}.{start}.{end}.csv"
    shell:
        """
        python3 workflow/scripts/query_gwas_catalog.py {wildcards.group} {wildcards.phenotype} {wildcards.chrom} {wildcards.start} {wildcards.end}
        """


# rule query_all_vep:
#     input:
#         chunk_list
#     output:
#         "output/meta-analysis/query/vep_all.csv"
#     run:
#         concat_list = [pd.read_csv(f) for f in input]
#         out = pd.concat(concat_list)
#         out.sort_values(['chrom', 'ps'], inplace = True)
#         out.to_csv(output[0], index = False)


# rule query_vep:
#     input: "output/meta-analysis/query/vep_chunks/chunk_{num}.txt"
#     output: "output/meta-analysis/query/vep_chunks/chunk_{num}.csv"
#     run:
#         data = json.dumps({'variants': pd.read_csv(input[0], header = None)[0].to_list()})
#         url = "https://rest.ensembl.org/vep/homo_sapiens/region"
#         headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
#         r = requests.post(url, headers = headers, data = data)

#         data = r.json()
#         rows = list()
#         for var in data:
#             chrom, pos, _, ref, alt, *_ = var['input'].split()
#             consq = var['most_severe_consequence']            
#             rows.append([chrom, pos, ref, alt, consq])
#         out = pd.DataFrame(rows, columns = ['chrom', 'ps', 'a1', 'a2', 'ensembl_consequence'])
#         out.to_csv(output[0], header = True, index = False)


"""
Snakefile 7. GWAS curation

"""
import yaml
import numpy as np
import pandas as pd
include: "read-config.smk"


peaklist = f"output/meta-analysis/peaks/peaklist/all.{config['group']}.peaklist"
peaklist = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
gwas_query_list = [f'output/meta-analysis/query/gwas/{r.group}.{r.phenotype}.{r.chrom}.{r.start}.{r.end}.csv' for _, r in peaklist.iterrows()]


rule assess_prev_assoc_signals:
    input:
        peaks=gwas_query_list,
        yaml="output/meta-analysis/query/gwas.filter.yaml"
    output:
        "output/meta-analysis/query/gwas.assessed.csv"
    run:
        concat_list = [pd.read_csv(peak) for peak in input.peaks]
        out = pd.concat(concat_list).reset_index(drop=True)

        with open(input.yaml, 'r') as f:
            d = yaml.safe_load(f)
            
        out['gwas previously associated'] = False

        for pheno, keywords in d.items():
            out.loc[(out['phenotype']==pheno)
                    & ((out['trait'].isin(keywords))
                        | (out['pval_desc_left'].isin(keywords))
                        | (out['pval_desc_right'].isin(keywords))
                        ),
                        'gwas previously associated'] = True

        prev_assocs = out.loc[out['gwas previously associated'], ['phenotype', 'rsid', 'gwas previously associated']].drop_duplicates()
        rs_pos = out[['rsid', 'chrom', 'start', 'end']].drop_duplicates()
        out = out[['phenotype', 'rsid', 'trait', 'pval_desc']].drop_duplicates()

        groupby = out.groupby(['phenotype', 'rsid'], as_index=False).agg(lambda x: ';'.join(x))
        groupby['trait'] = [';'.join(set(i.split(';'))) for i in groupby['trait']]
        groupby['pval_desc'] = [';'.join(set(i.split(';'))) for i in groupby['pval_desc']]
        groupby = groupby.merge(rs_pos, how = 'left')
        groupby = groupby.merge(prev_assocs, how = 'left')
        groupby['gwas previously associated'].fillna(False, inplace = True)

        groupby.rename(columns = {
                                'rsid': 'ensembl_rs',
                                'trait': 'gwas associated trait',
                                'pval_desc': 'gwas associated trait description'
                                }, inplace = True)
        groupby.drop_duplicates(inplace = True)
        groupby.to_csv(output[0], index = False)
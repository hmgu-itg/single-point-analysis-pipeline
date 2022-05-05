import re
import sys
import json

import requests
import pandas as pd


def snp_search(chrom, start, end):
    endpoint = 'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByChromBpLocationRange'
    params = f'?chrom={chrom}&bpStart={start}&bpEnd={end}'
    params2 = '&page={}&size=500'.format

    url = endpoint + params
    r = requests.get(url + params2(0))
    d = r.json()
    page_info = d['page']
    
    queries = [d]
    if page_info['totalPages'] > 1:
        print(f'Multiple pages for {chrom}:{start}-{end}')
        

        for i in range(1, page_info['totalPages']):
            r = requests.get(url + params2(i))
            d = r.json()
            page_info = d['page']
            queries.append(d)

    q_data = [var for q in queries for var in q['_embedded']['singleNucleotidePolymorphisms']]
    return q_data


def snp_associations(rsid: str):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations?projection=associationBySnp"
    r = requests.get(url)
    d = r.json()
    rows = list()
    for assoc in d['_embedded']['associations']:
        pval = assoc['pvalue']
        trait = ';'.join([t['trait']for t in assoc['efoTraits']])
        pval_desc = assoc['pvalueDescription']
        
        rows.append([rsid, pval, trait, pval_desc])
    return pd.DataFrame(rows, columns = ['rsid', 'pval', 'trait', 'pval_desc'])


def process_associations(associations):
    new = associations.copy()
    # new['pval_desc'] = [re.search(r'(?<=\().*(?=\))', i).group() if isinstance(i, str) else None for i in new['pval_desc']]
    new_pval_desc = list()
    for i in new['pval_desc']:
        if isinstance(i, str):
            m = re.search(r'(?<=\().*(?=\))', i)
            if m is not None:
                new_pval_desc.append(m.group())
            else:
                i = i.strip()
                if i == '':
                    new_pval_desc.append(None)
                else:
                    new_pval_desc.append(i)
        else:
            new_pval_desc.append(None)
    new['pval_desc'] = new_pval_desc
    
    
    rows = list()
    for i in new['pval_desc']:
        if isinstance(i, str):
            if 'mitochondrial' in i and i.count(',')>1:
                *left, right = i.split(',')
                left = ','.join(left).strip()
                right = right.strip()

                rows.append([i, left, right])
            elif i.count(',')==1:
                left, right = i.split(',')
                rows.append([i, left.strip(), right.strip()])

            elif i is None:
                rows.append([i, None, None])

            else:
                rows.append([i, i, None])
    
    
    left_right = pd.DataFrame(rows, columns = ['pval_desc', 'pval_desc_left', 'pval_desc_right']).drop_duplicates()

    new = new.merge(left_right, on = 'pval_desc')
    return new.reset_index(drop = True)


def process_pval_desc_right(series):
    r = list()
    for i in series:
        if isinstance(i, str) and re.match(r'.*(\.\d+)+$', i) is not None:
            r.append(re.sub(r'(\.\d+)+$', '', i))
        else:
            r.append(i)
    return r


def get_variant_positions(rs_list: list):
    data_dict = dict()
    maxlimit = 200
    for i in range(0, len(rs_list), maxlimit):
        data = json.dumps({'ids': rs_list[i:i+maxlimit]})
        url = "https://rest.ensembl.org/variation/homo_sapiens"
        headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
        r = requests.post(url, headers=headers, data=data)

        if not r.ok:
            r.raise_for_status()
        data_dict.update(r.json())

    rows = list()
    for rs, info in data_dict.items():
        mappings = info['mappings']
        if len(mappings)>1:
            for m in mappings:
                if m['assembly_name']=='GRCh38':
                    mapping = m
                    break
            else:
                print(f'[WARNING] No GRCh38 mapping info for {rs}!')
                mapping = mappings[0]
                print(f'[WARNING] Using assembly mapping from {mapping["assembly_name"]}')
        else:
            mapping = mappings[0]
        chrom = mapping['seq_region_name']
        start = mapping['start']
        end = mapping['end']
        rows.append((rs, chrom, start, end))
    
    return pd.DataFrame(rows, columns = ['rsid', 'chrom', 'start', 'end'])


def get_associations(chrom, start, end):
    q_data = snp_search(chrom, start, end)

    rs_ids = {re.sub(r'\/.*', '', q['rsId']) for q in q_data}
    snp_assocs = list()
    for rs in rs_ids:
        try:
            snp_assoc = snp_associations(rs)
        except json.JSONDecodeError:
            print(f'[WARNING] Failed to get associations for {rs}')
            continue
        snp_assocs.append(snp_assoc)
    
    try:
        associations = pd.concat(snp_assocs)
    except ValueError:
        return pd.DataFrame(columns = ['rsid', 'chrom', 'start', 'end', 'pval', 'trait', 'pval_desc', 'pval_desc_left', 'pval_desc_right'])
    associations = process_associations(associations)
    associations['pval_desc_right'] = process_pval_desc_right(associations['pval_desc_right'])
    var_positions = get_variant_positions(associations['rsid'].to_list())
    associations = associations.merge(var_positions, on = 'rsid', how = 'left')
    associations = associations[['rsid', 'chrom', 'start', 'end', 'pval', 'trait', 'pval_desc', 'pval_desc_left', 'pval_desc_right']]
    return associations



if __name__ == '__main__':
    _, group, phenotype, chrom, start, end = sys.argv

    df = get_associations(chrom, start, end)

    df.insert(0, 'phenotype', phenotype)
    df.insert(0, 'group', group)
    df.to_csv(f"output/meta-analysis/query/gwas/{group}.{phenotype}.{chrom}.{start}.{end}.csv", index = False, doublequote=True)
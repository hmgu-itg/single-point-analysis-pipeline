"""
Snakefile 4. Plot peaks

Run PeakPlotter for all detected peaks

Example
-------
# To run peakplotter on all detected peaks
$ snakemake --cores 10 --restart-times 3 --keep-going --use-singularity --snakefile workflow/rules/peakplot.smk 

# To run a single peakplotter run
$ snakemake --cores 1 --restart-times 3 --keep-going --use-singularity --snakefile workflow/rules/peakplot.smk output/meta-analysis/peaks/{group}.{phenotype}/{chrom}.{start}.{end}.500kb.csv
"""
from pathlib import Path

import pandas as pd
include: "read-config.smk"
include: "meta-analysis.smk"


rule get_vep_query_list:
    input:
        "output/meta-analysis/peaks/all.peaks.filtered.csv.gz"
    output:
        "output/meta-analysis/query/vep_chunks.txt",
        directory("output/meta-analysis/query/vep_chunks")
    run:
        peaks = pd.read_csv(input[0])
        novel = peaks['ensembl_consequence']=='novel'
        subset = peaks.loc[novel, ['chrom', 'ps', 'a1', 'a2']].drop_duplicates()

        # Create POST input data
        variants = [f"{chrom} {pos} . {ref}  {alt}  . . ." for _, (chrom, pos, ref, alt) in subset.iterrows()]
        max_limit = 200
        chunks = [variants[i:i+max_limit] for i in range(0, len(variants), max_limit)]
        path = Path(output[1])
        path.mkdir(parents=True, exist_ok=True)
        files = list()
        for idx, chunk in enumerate(chunks):
            filename = path.joinpath(f"chunk_{idx}.txt")
            files.append(filename)
            with open(filename, 'w') as f:
                for line in chunk:
                    f.write(line)
                    f.write('\n')
        pd.DataFrame(files).to_csv(output[0], index = False, header = False)


def collect_all_peak_csvs_input(w):
    peaklist = f"output/meta-analysis/peaks/peaklist/all.{config['group']}.peaklist"

    all_peaks = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
    peaks = [f'output/meta-analysis/peaks/{row.group}.{row.phenotype}/{row.chrom}.{row.start}.{row.end}.500kb.csv' for _, row in all_peaks.iterrows()]
    return peaks

rule collect_all_peak_csvs:
    input:
        collect_all_peak_csvs_input
    params:
        p_threshold=config['QC_thresholds']['p-value']
    output:
        "output/meta-analysis/peaks/all.peaks.csv.gz",
        "output/meta-analysis/peaks/all.peaks.filtered.csv.gz"
    run:
        concat_list = list()
        for f in input:
            string = re.sub(r'.*\/peaks\/', '', f)
            string = re.sub(r'\/.*', '', string)
            group, phenotype = string.split('.')
            peak = f.split('/')[-1].replace('.500kb.csv', '')
            df = pd.read_csv(f, header = 0)
            df.insert(0, 'peak', peak)
            df.insert(0, 'phenotype', phenotype)
            df.insert(0, 'group', group)
            concat_list.append(df)
        all_data = pd.concat(concat_list)
        all_data.to_csv(output[0], index = False, header = True, compression='gzip')
        all_data[all_data['p-value']<=params.p_threshold].to_csv(output[1], index = False, header = True, compression='gzip')


rule peak_metal:
    input:
        # metal="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz",
        # bfiles=multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam')
        metal=rules.filter_metal.output.gz,
        bfiles=rules.phenotype_mac_filter.output.bfile
    params:
        # bfile="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        bfile=rules.phenotype_mac_filter.params.out,
        outdir="output/meta-analysis/peaks/{group}.{phenotype}",
        chrom="{chrom}",
        start="{start}",
        end="{end}"
    singularity: "library://hmgu-itg/default/peakplotter:0.4.3"
    output:
        "output/meta-analysis/peaks/{group}.{phenotype}/{group}.{phenotype}.{chrom}.{start}.{end}.log",
        multiext("output/meta-analysis/peaks/{group}.{phenotype}/{chrom}.{start}.{end}.500kb", '.html', '.csv')
    resources:
        mem_per_cpu=lambda w: "5G" if (int(w.end) - int(w.start)) > 5_000_000 else "2G"
    shell:
        """
        peakplotter-region  \
          --chr-col Chrom \
          --pos-col Pos \
          --rs-col MarkerName \
          --pval-col P-value \
          --a1-col Allele1 \
          --a2-col Allele2 \
          --maf-col Freq1 \
          --build 38 \
          --assoc-file {input.metal} \
          --bfiles {params.bfile} \
          --out {params.outdir} \
          --chrom {params.chrom} \
          --start {params.start} \
          --end {params.end} \
          --debug
        """

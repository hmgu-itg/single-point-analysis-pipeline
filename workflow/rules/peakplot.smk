"""
Snakefile 4. Plot peaks

Run PeakPlotter for all detected peaks

Example
-------
# To run peakplotter on all detected peaks
$ snakemake --cores 10 --restart-times 3 --keep-going --use-singularity --snakefile workflow/rules/peakplot.smk 

# To run a single peakplotter run
$ snakemake --cores 1 --restart-times 3 --keep-going --use-singularity --snakefile workflow/rules/peakplot.smk output/meta-analysis/peaks/{group}.{phenotype}/{chrom}.{start}.{end}.500kb
"""
from pathlib import Path

import pandas as pd
include: "read-config.smk"

peaklist = expand("output/meta-analysis/peaks/peaklist/all.{group}.peaklist", group=config['group'])[0]

all_peaks = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
peaks = [f'{row.group}.{row.phenotype}/{row.chrom}.{row.start}.{row.end}' for _, row in all_peaks.iterrows()]

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

rule collect_all_peak_csvs:
    input:
        expand("output/meta-analysis/peaks/{peak}.500kb.csv", peak = peaks)
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
        metal="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz",
        bfiles=multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam')
    params:
        bfile="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        outdir="output/meta-analysis/peaks/{group}.{phenotype}",
        chrom="{chrom}",
        start="{start}",
        end="{end}"
    singularity: "library://hmgu-itg/default/peakplotter:dev"
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
          --a1-col Ref \
          --a2-col Alt \
          --maf-col Alt_Freq \
          --build 38 \
          --assoc-file {input.metal} \
          --bfiles {params.bfile} \
          --out {params.outdir} \
          --chrom {params.chrom} \
          --start {params.start} \
          --end {params.end} \
          --debug
        """

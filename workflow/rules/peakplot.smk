import pandas as pd
include: "read-config.smk"

peaklist = expand("output/meta-analysis/peaks/peaklist/all.{group}.peaklist", group=config['group'])[0]

all_peaks = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
peaks = [f'{row.group}.{row.phenotype}/{row.chrom}.{row.start}.{row.end}' for _, row in all_peaks.iterrows()]

rule all4:
    input:
        expand("output/meta-analysis/peaks/{peak}.500kb.html", 
                peak = peaks)

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
    singularity: "library://hmgu-itg/default/peakplotter"
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
          --end {params.end}
        """

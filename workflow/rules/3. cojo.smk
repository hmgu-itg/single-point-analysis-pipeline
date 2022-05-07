"""
Snakefile 5. Run COJO

Example
-------
# To run cojo on all detected peaks
$ snakemake --cores 10 --keep-going --use-singularity --snakefile workflow/rules/cojo.smk run_all_cojo

"""
include: "2. single-cohort.smk"

def run_all_cojo_input(w):
    peaklist = checkpoints.detect_peaks.get(cohort=w.cohort, group=w.group, phenotype=w.phenotype).output[0]
    try:
        peaklist = pd.read_csv(peaklist, sep = '\t', header = None, names = ['group', 'phenotype', 'chrom', 'start', 'end'])
    except pd.errors.EmptyDataError:
        return []
    peaks = [f'output/{w.cohort}/{group}/{phenotype}/cojo/{chrom}.{start}.{end}.jma.cojo' for _, (group, phenotype, chrom, start, end) in peaklist.iterrows()]
    return peaks

rule run_all_cojo:
    input:
        run_all_cojo_input
    output:
        "output/{cohort}/{group}/{phenotype}/all.jma.csv.gz"
    run:
        if len(input)==0: # When no signif signal
            empty = pd.DataFrame(columns = ["group", "phenotype", "peak", "Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r", "cojo_tophit"])
            empty.to_csv(output[0], index = False, header = True, compression = 'gzip')
            return

        concat_list = list()
        for f in input:
            string = re.sub(r'.*\/cojo\/', '', f)
            string = re.sub(r'\/.*', '', string)
            peak = f.split('/')[-1].replace(f'{string}.', '').replace('.jma.cojo', '')
            df = pd.read_csv(f, sep = '\t', header = 0)
            min_p = df['p'].min()
            df['is cojo tophit'] = False
            df.loc[df['p'] == min_p, 'cojo_tophit'] = True
            df.insert(0, 'group', wildcards.group)
            df.insert(1, 'phenotype', wildcards.phenotype)
            df.insert(2, 'peak', peak)
            concat_list.append(df)
        pd.concat(concat_list).to_csv(output[0], index = False, header = True, compression = 'gzip')


rule make_cojofile:
    '''
    Creates cojofile from bfile and assoc file.
    bfile is only used to calculate the "sample size" value for each variant (The "N" column of cojofile).
    The input assoc file makes up bulk of the cojofile. 
    
    Reference
    --------
    1. https://yanglab.westlake.edu.cn/software/gcta/#COJO
    2. https://yanglab.westlake.edu.cn/software/gcta/#MLMA
    '''
    input:
        assoc=rules.gcta.output.mlma_bgz,
        samplesize=rules.get_samplesize.output.samplesize
    params:
        bfile=rules.filter_bfile.params.out,
        prefix="output/{cohort}/{group}/{phenotype}/cojofile/cojofile"
    output:
        "output/{cohort}/{group}/{phenotype}/cojofile/cojofile.ma.gz"
    run:
        samplesize = pd.read_csv(input.samplesize,
                                sep = ' ',
                                header = None,
                                names = ['SNP', 'N'])
        assoc = pd.read_csv(input.assoc,
                            sep = '\t',
                            usecols = ['SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p'])
        cojofile = assoc.merge(samplesize)
        assert assoc.shape[0] == cojofile.shape[0], 'association file and cojofile row count do not match!'
        cojofile.rename(columns = {'Freq': 'freq'}, inplace=True)

        cojofile = cojofile[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']] # Just to be absolutely sure
        cojofile.to_csv(output[0],
                        sep = ' ',
                        index = False,
                        compression = 'gzip')



rule cojo:
    input:
        bfiles=rules.filter_bfile.output,
        cojofile=rules.make_cojofile.output,
        assoc=rules.gcta.output.mlma_bgz
    params:
        bfile=rules.filter_bfile.params.out,
        threshold=config['QC_thresholds']['p-value'],
        prefix="output/{cohort}/{group}/{phenotype}/cojo/{chrom}.{start}.{end}"
    output:
        multiext("output/{cohort}/{group}/{phenotype}/cojo/{chrom}.{start}.{end}", 
            ".jma.cojo",
            ".cma.cojo",
            ".ldr.cojo")
    log:
        "output/{cohort}/{group}/{phenotype}/cojo/{chrom}.{start}.{end}.cojo.log"
    shell:
        """
        workflow/scripts/run_cojo.sh \
            {params.bfile} \
            {input.cojofile} \
            {input.assoc} \
            {params.threshold} \
            {wildcards.chrom} \
            {wildcards.start} \
            {wildcards.end} \
            {params.prefix} \
            {resources.mem_mb} 2>&1 > {log}
        """

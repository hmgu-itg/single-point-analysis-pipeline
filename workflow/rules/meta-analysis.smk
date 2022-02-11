"""
Snakefile 2.

$ snakemake --cores 100 --snakefile workflow/rules/meta-analysis.smk --use-singularity --batch create_all_metal=1/10
$ snakemake --cores 100 --snakefile workflow/rules/meta-analysis.smk --use-singularity create_all_metal
"""
include: "read-config.smk"
container: config['container']


# TODO: How should we deal with variants that's only in one cohort?

rule create_all_metal:
    input:
        expand("output/meta-analysis/metal/{group}.{phenotype}.filtered.txt.gz", group = config['group'], phenotype = config['phenotypes']),
        "output/bfile/combined.frq2"

rule metal:
    """
    METAL analysis wrapper for GCTA files as input
    
    Note
    ----
    A1 is the effect allele and A2 is the other allele
    in GCTA output files, and METAL documentation suggest
    putting the effect allele first.
    """
    version: '0.2'
    input:
        lambda w: expand("output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz", \
               cohort = config['cohorts'], allow_missing=True)
    params:
        prefix="output/meta-analysis/metal/{group}.{phenotype}"
    output:
        metal="output/meta-analysis/metal/{group}.{phenotype}.txt.gz",
        index="output/meta-analysis/metal/{group}.{phenotype}.txt.gz.tbi",
        info="output/meta-analysis/metal/{group}.{phenotype}.info",
        cmd=temp("output/meta-analysis/metal/{group}.{phenotype}.cmd")
    log:
        "output/meta-analysis/metal/{group}.{phenotype}.log"
    resources:
        cpus_per_task=1,
        mem_mb=2000,
        time="1:0:0",
    shell:
        """
        # Create command file
        echo 'SEPARATOR TAB'  > {output.cmd}
        echo 'MARKER SNP'     >> {output.cmd}
        echo 'ALLELE A1 A2'   >> {output.cmd} 
        echo 'FREQLABEL Freq' >> {output.cmd}
        echo 'EFFECT b'       >> {output.cmd}
        echo 'STDERR se'      >> {output.cmd}
        echo 'PVALUE p'       >> {output.cmd}
        echo 'SCHEME STDERR'  >> {output.cmd}
        echo 'AVERAGEFREQ ON' >> {output.cmd}
        echo 'MINMAXFREQ ON'  >> {output.cmd}
        for f in {input} ; do
          echo \"PROCESS $f\" >> {output.cmd}
        done
        echo 'OUTFILE {params.prefix}. .txt' >> {output.cmd}
        echo 'ANALYZE' >> {output.cmd}

        metal < {output.cmd} 2>&1 > {log}

        ## CLEAN, COMPRESS, INDEX
        cat \
          <(echo "Chrom MarkerName Pos Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr P-value Direction" | tr ' ' '\\t') \
          <(tail -n+2 {params.prefix}.1.txt | sed 's/^chr//' | tr ':' '\\t' | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,\"chr\"$1\":\"$2,$2,toupper($3),toupper($4),$5,$6,$7,$8,$9,$10,$11,$12}}' | sort -k1,1n -k3,3n) \
          | bgzip > {output.metal} 2>> {log}

        tabix -s1 -b3 -e3 -S1 {output.metal} 2>> {log}
        rm {params.prefix}.1.txt

        ## Modify info file
        mv {params.prefix}.1.txt.info {output.info}

        sed -i 's,{params.prefix}.1.txt,{output.metal},' {output.info}
        sed -i 's/Marker    /MarkerName  /' {output.info}
        """



rule merge_bfiles:
    """
    Merges all cohort genotype bfiles into a single bfile.

    Outputs
    -------
    bfile: All cohort bfiles combined into one bfile
    missnp: A list of variant IDs of multiallelic variants excluded from the combined bfile
    """
    input:
        expand("output/bfile/{cohort}.{ext}", cohort=config['cohorts'], ext=['bed', 'bim', 'fam'])
    params:
        input=expand("output/bfile/{cohort}", cohort=config['cohorts']),
        out="output/bfile/combined"
    output:
        bfile=multiext("output/bfile/combined", '.bed', '.bim', '.fam'),
        missnp="output/bfile/combined-merge.missnp"
    threads: 20
    resources:
        cpus_per_task=20,
        mem_mb=5000
    shell:
        """
        mergelist=$(mktemp /tmp/single-point-analysis-pipeline.merge_bfiles.XXXXXX)
        exec 3>\"$mergelist\"
        echo {params.input} | tr ' ' '\\n' >&3

        plink --allow-no-sex --threads {threads} --make-bed \
          --merge-list $mergelist \
          --out {params.out} || true # || true to prevent exitting due to bash strict mode

        # If merge failed due to multiallelics, exclude multiallelics and re-merge
        if [ -f "{output.missnp}" ]
        then
          exec 3>\"$mergelist\"
          for bfile in {params.input}
          do
            tmp_bfile=${{bfile}}-combine-tmp
            plink --allow-no-sex --threads {threads} --make-bed \
              --bfile $bfile \
              --exclude {output.missnp} \
              --out $tmp_bfile
            echo $tmp_bfile >&3
          done
          
          plink --allow-no-sex --threads {threads} --make-bed \
            --merge-list $mergelist \
            --out {params.out}
        else
          touch {output.missnp}
        fi
        rm output/bfile/*-combine-tmp.*
        exec 3>-
        """

rule combined_freq:
    input:
        multiext("output/bfile/combined", '.bed', '.bim', '.fam')
    params:
        bfile="output/bfile/combined",
        awk="{if (NR!=1){print $3,$6}}" # Variant ID and MAF column
    output:
        frq="output/bfile/combined.frq",
        frq2="output/bfile/combined.frq2"
    threads: 20
    resources:
        cpus_per_task=20,
        mem_per_cpu="5G"
    shell:
        """
        plink --threads {threads} \
            --bfile {params.bfile} \
            --freq \
            --out {params.bfile}
        
        awk -F '[[:space:]]+' '{params.awk}' {output.frq} > {output.frq2}
        """

rule phenotype_mac_filter:
    """
    `output/bfile/mac/all.mac.txt` file is created in the `variant-qc.smk` snakefile.
    
    MAC==10 equivalent of MAF is calculated by now, and this is the filtering rule. 
    1. The MAC==10 equivalent MAF for the phenotype is extracted from the `all.mac.txt` file
    2. Then a `{params.out}.mac.excludelist` file is created which contains all variant IDs with 10 or less MAC
    3. The excludelist file is then used as input to plink to remove those variants from the bfile.

    Output
    ------
    bfile: 
    excludelist: List of variant IDs with less than the MAC threshold
    """
    input:
        mac="output/bfile/mac/all.mac.txt",
        bfiles=multiext("output/bfile/combined", '.bed', '.bim', '.fam', '.frq2')
    params:
        bfile="output/bfile/combined",
        out="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        id="{group}.{phenotype}"
    output:
        bfile=multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam', '.nosex'),
        excludelist="output/meta-analysis/temp-bfiles/{group}.{phenotype}.mac.excludelist"
    threads: 20
    resources:
        cpus_per_task=20,
        mem_mb=5000,
        time="2:0:0"
    shell:
        """
        mac_threshold=$(awk -F' ' -v id='{params.id}' '{{if ($1==id){{print $4}}}}' {input.mac})

        awk -F ' ' -v mac_threshold=\"$mac_threshold\" '{{if ($2<=mac_threshold){{print $1}}}}' {params.bfile}.frq2 > {output.excludelist}

        plink --allow-no-sex --threads {threads} --make-bed \
              --bfile {params.bfile} \
              --exclude {output.excludelist} \
              --out {params.out}
        """


rule filter_metal:
    input:
        metal=rules.metal.output.metal,
        missnp=rules.merge_bfiles.output.missnp,
        excludelist=rules.phenotype_mac_filter.output.excludelist
    output:
        gz="output/meta-analysis/metal/{group}.{phenotype}.filtered.txt.gz",
        tbi="output/meta-analysis/metal/{group}.{phenotype}.filtered.txt.gz.tbi"
    shell:
        """
        zgrep -v -w -f <(cat {input.missnp} {input.excludelist}) {input.metal} | grep -v na | bgzip > {output.gz}
        tabix -s1 -b3 -e3 -S1 {output.gz}
        """


rule manqq_metal:
    input: rules.metal.output.metal # "output/meta-analysis/metal/{group}.{phenotype}.metal.gz"
    params:
        out_prefix="output/meta-analysis/manqq/{group}.{phenotype}",
        maf=0.0
    output:
        multiext("manqq/{group}.{phenotype}", 
                 '.run_conf',
                 '.qq.png',
                 '.lambda.txt',
                 '.man.png')
    singularity: "library://hmgu-itg/default/manqq"
    log:
        "output/meta-analysis/manqq/{group}.{phenotype}.log"
    shell:
        """
        run_manqq.R \
        --chr-col Chrom \
        --pval-col P-value \
        --pos-col Pos \
        --a1 Alt \
        --a2 Ref \
        --build 38 \
        --image png \
        --af-col Alt_Freq \
        --maf-filter {params.maf}
        {input} \
        {params.out_prefix} 2>&1 > {log}
        """

use rule manqq_metal as manqq_metal_filtered with:
    params:
        out_prefix="output/meta-analysis/manqq/{group}.{phenotype}.filtered",
        maf=0.001
    output:
        multiext("manqq/{group}.{phenotype}.filtered", 
                 '.run_conf',
                 '.qq.png',
                 '.lambda.txt',
                 '.man.png')
    log:
        "output/meta-analysis/manqq/{group}.{phenotype}.filtered.log"
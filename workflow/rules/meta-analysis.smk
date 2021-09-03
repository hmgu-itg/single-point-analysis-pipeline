include: "read-config.smk"
container: config['container']

rule all3:
    input:
        expand("output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz", \
               cohort = config['cohorts'], group = config['group'], phenotype = config['phenotypes'])


rule run_metal:
    input:
        lambda w: expand("output/single-cohort/gcta/{cohort}/{cohort}.{group}.{phenotype}.mlma.gz", \
               cohort = config['cohorts'], allow_missing=True)
    params:
        out_prefix="output/meta-analysis/metal/{group}.{phenotype}"
    output:
        ver="output/meta-analysis/metal/{group}.{phenotype}.metal.ver",
        info="output/meta-analysis/metal/{group}.{phenotype}.metal.info",
        bgz="output/meta-analysis/metal/{group}.{phenotype}.metal.gz",
        tbi="output/meta-analysis/metal/{group}.{phenotype}.metal.gz.tbi"
    log:
        "output/meta-analysis/metal/{group}.{phenotype}.metal.log"
    shell:
        """
        ## METAL analysis wrapper for GCTA files as input
        metal_cmd={params.out_prefix}.metal.cmd
        metal --version | head -2 > {params.out_prefix}.metal.ver

        # WARNING: Notice that the order of allele is A2 then A1?
        # This is because A2 is the reference allele and A1 is the effect allele in
        # GCTA output files.
        ## RUN METAL
        echo 'SEPARATOR TAB'  > $metal_cmd
        echo 'MARKER SNP'     >> $metal_cmd
        echo 'ALLELE A2 A1'   >> $metal_cmd
        echo 'FREQLABEL Freq' >> $metal_cmd
        echo 'EFFECT b'       >> $metal_cmd
        echo 'STDERR se'      >> $metal_cmd
        echo 'PVALUE p'       >> $metal_cmd
        echo 'SCHEME STDERR'  >> $metal_cmd
        echo 'AVERAGEFREQ ON' >> $metal_cmd
        echo 'MINMAXFREQ ON'  >> $metal_cmd
        for f in {input} ; do
          echo \"PROCESS $f\" >> $metal_cmd
        done
        echo 'OUTFILE {params.out_prefix}. .txt' >> $metal_cmd
        echo 'ANALYZE' >> $metal_cmd

        metal < $metal_cmd 2>&1 > {log}

        ## CLEAN, COMPRESS, INDEX
        cat \
        <(echo "Chrom MarkerName Pos Ref Alt Alt_Freq FreqSE MinFreq MaxFreq Effect StdErr P-value Direction" | tr ' ' '\\t') \
        <(tail -n+2 {params.out_prefix}.1.txt | sed 's/^chr//' | tr ':' '\\t' | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,\"chr\"$1\":\"$2,$2,toupper($4),toupper($3),$5,$6,$7,$8,$9,$10,$11,$12}}' | sort -k1,1n -k3,3n) \
        | bgzip > {output.bgz} 2>> {log}
        tabix -s1 -b3 -e3 -S1 {output.bgz} 2>> {log}
        sed 's,{params.out_prefix}.1.txt,{output.bgz},' {params.out_prefix}.1.txt.info |
          sed 's/Marker    /MarkerName  /' |
          sed 's/Allele1/Alt    /' |
          sed 's/Allele2/Ref    /' |
          sed 's/Freq1     /Alt_Freq/' |
          sed 's/allele *1/Alt allele/' > {output.info}
        rm {params.out_prefix}*.txt* $metal_cmd 2>> {log}
        """


rule merge_bfiles:
    input:
        expand("output/bfile/{cohort}.{ext}", cohort=config['cohorts'], ext=['bed', 'bim', 'fam'])
    params:
        input=expand("output/bfile/{cohort}", cohort=config['cohorts']),
        out="output/bfile/combined"
    output:
        multiext("output/bfile/combined", '.bed', '.bim', '.fam'),
        missnp="output/bfile/combined-merge.missnp"
    threads: 20
    resources:
        cpus_per_task=20,
        mem_per_cpu="5G"
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
        awk="{if (NR!=1){print $3,$6}}"
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
    input:
        mac="output/bfile/mac/all.mac.txt",
        bfiles=multiext("output/bfile/combined", '.bed', '.bim', '.fam', '.frq2')
    params:
        bfile="output/bfile/combined",
        out="output/meta-analysis/temp-bfiles/{group}.{phenotype}",
        id="{group}.{phenotype}"
    output:
        multiext("output/meta-analysis/temp-bfiles/{group}.{phenotype}", '.bed', '.bim', '.fam', '.nosex', '.mac.excludelist')
    threads: 20
    resources:
        cpus_per_task=20,
        mem_per_cpu="5G"
    shell:
        """
        mac_threshold=$(awk -F' ' -v id='{params.id}' '{{if ($1==id){{print $4}}}}' {input.mac})
        excludelist={params.out}.mac.excludelist

        awk -F ' ' -v mac_threshold=\"$mac_threshold\" '{{if ($2<=mac_threshold){{print $1}}}}' {params.bfile}.frq2 > $excludelist

        plink --allow-no-sex --threads {threads} --make-bed \
              --bfile {params.bfile} \
              --exclude $excludelist \
              --out {params.out}
        """

rule filter_metal:
    input:
        metal="output/meta-analysis/metal/{group}.{phenotype}.metal.gz",
        missnp="output/bfile/combined-merge.missnp",
        excludelist="output/meta-analysis/temp-bfiles/{group}.{phenotype}.mac.excludelist"
    output:
        gz="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz",
        tbi="output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz.tbi"
    shell:
        """
        zgrep -v -w -f <(cat {input.missnp} {input.excludelist}) {input.metal} | grep -v na | bgzip > {output.gz}
        tabix -s1 -b3 -e3 -S1 {output.gz}
        """


rule collect_peaks:
    input:
        "output/meta-analysis/temp-bfiles/{group}.{phenotype}.metal.filtered.gz"
    params:
        span=config['peakplotter']['span'],
        signif=5e-8,
        group="{group}",
        phenotype="{phenotype}"
    output:
        "output/meta-analysis/peaks/peaklist/{group}.{phenotype}.peaklist"
    singularity: "library://hmgu-itg/default/peakplotter"
    shell:
        "python3 workflow/scripts/collect_peaks.py {input} {params.span} {params.signif} {params.group} {params.phenotype} {output}"
        

rule collect_all_peaks:
    input:
        expand("output/meta-analysis/peaks/peaklist/{group}.{phenotype}.peaklist", 
                phenotype = config['phenotypes'], allow_missing=True)
    output:
        "output/meta-analysis/peaks/peaklist/all.{group}.peaklist"
    run:
        all_df = list()
        for f in input:
            try:
                df = pd.read_csv(f, sep = '\t', header = None)
            except pd.errors.EmptyDataError:
                continue
            all_df.append(df)
        all_df = pd.concat(all_df)
        all_df.sort_values([0, 1, 2, 3, 4], inplace = True)
        all_df.to_csv(output[0], sep = '\t', header = False, index = False)


rule manqq_metal:
    input: "output/meta-analysis/metal/{group}.{phenotype}.metal.gz"
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
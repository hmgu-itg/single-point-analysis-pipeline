# Prerequisites
1. Snakemake
2. Singularity

# Input data
- Plink binary genotype data (`.bed`, `.bim`, `.fam`)
    - See https://www.cog-genomics.org/plink/1.9/input#bed
    - Recommended to filter out variants with extremely high missingness[^1].
- Plink `.lmiss` file
    - Created from the binary genotype data
- GCTA-GRM data (`.grm.bin`, `.grm.id`, `.grm.N.bin`)
    - Created from the binary genotype data
- Phenotype values (text file)

<details>
<summary>More information about input data</summary>

### Plink `.lmiss` file  
See https://www.cog-genomics.org/plink/1.9/formats#lmiss

```bash
plink \
    --bfile input-bfile \
    --memory 1000 \
    --threads 1 \
    --missing \
    --out output-bfile \
```
Please refer to the [System resource usage](https://www.cog-genomics.org/plink/1.9/other#memory) section of Plink for more information about `--memory` and `--threads`.

### GCTA-GRM data
See https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM

### Phenotype values
The format of the phenotype file needs to be as shown below. 

|||
|---|---|
|SampleA|0.312|
|SampleB|0.150|

The file needs to be a headerless, tab-separated file with two columns: first column with sample ID that matches the sample IDs in the `.fam` file, second column with the phenotype value for the corresponding samples.
</details>

# Running the analysis
## Before running
Examine the following two configuration files:
- `config.yaml`
- `analyse/config.yaml`

The `config.yaml` file is used to configure the can be modified to is where the user can define configuration values which will affect the overall result of the analysis. The `analyse.config.yaml` is a [snakemake configuration profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for running the workflow. Modify the values here based on your server resource availability.

## Running the pipeline
To run the pipeline, you must specify the `bfile`, `grm`, `lmiss`, and `phenotype_file` configuration values either in the `config.yaml` file or using the `--config` option (See [Snakemake Configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)).

For example, running one association analysis using the `--config` option:
```bash
snakemake \
    --profile analyse \
    --config \
        bfile=/path/to/cohort-bfile \
        lmiss=/path/to/cohort-bfile.lmiss \
        grm=/path/to/GCTA-grm \
        phenotype_file=/path/to/phenotype_file.txt \
    output/done
```

## Running the pipeline for multiple phenotypes
You can run the above command in a loop with different `phenotype` and `phenotype_file` value to run association on multiple traits.  
Alternatively, if you have a job scheduler available (e.g. Slurm), submit each individual run as jobs.  
If you are running multiple jobs in parallel, you may need to add the `--nolock` option to your snakemake command. The output files should not overlap as long as you specify a different combination of the `group` and `phenotype` values, but please be aware of [why Snakemake locks the working directory](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-lock-the-working-directory).

-----  

[^1]: This is due to the GCTA-COJO algorithm having issues analysing variants with low sample size.
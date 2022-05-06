# Prerequisites
1. Snakemake
2. Singularity

# Input data
- Plink binary genotype data (`.bed`, `.bim`, `.fam`)
    - See https://www.cog-genomics.org/plink/1.9/input#bed
- GCTA-GRM data (`.grm.bin`, `.grm.id`, `.grm.N.bin`)
    - See https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
- Phenotype values

## Phenotype values
The format of the phenotype file needs to be as shown below. 

|||
|---|---|
|SampleA|0.312|
|SampleB|0.150|

The file needs to be a tab-separated file with two columns: first column with sample ID that matches a sample ID in the `.fam` file, second column with the phenotype value for the corresponding sample.

# Running the analysis
## Before running
Examine the following two configuration files:
- `config.yaml`
- `analyse/config.yaml`

The first file contain configurations such as QC values which will affect the overall result of the analysis. The second file contain configurations for running the workflow. Modify the values here based on your server resource availability.

## Running one analysis
If you're only running association analysis for one trait, use the command below:

```bash
snakemake \
    run \
    --profile analyse \
    --config \
        cohort=MANOLIS \
        bfile=/path/to/cohort-bfile \
        grm=/path/to/GCTA-grm \
        group=Anthropometric \
        phenotype=BMI \
        phenotype_file=/path/to/phenotype_file.txt
```

## Running multiple analysis for a single cohort
If you want to run association analysis for multiple traits for a single cohort,  
you should first process the genotype file which will be shared across all association runs.
```bash
snakemake \
    filter_genotype \
    --profile analyse 
    --config \
        cohort=MANOLIS \
        bfile=/path/to/cohort-bfile
```
Once the above completes, run each individual association analysis:
```bash
snakemake \
    run \
    --profile analyse \
    --config \
        cohort=MANOLIS \
        grm=/path/to/GCTA-grm \
        group=Anthropometric \
        phenotype=BMI \
        phenotype_file=/path/to/phenotype_file.txt
```
You can run the above command in a loop with different `phenotype` and `phenotype_file` value to run association on all traits.  
Alternatively, submit each individual run as jobs to cluster to run all in parallel.

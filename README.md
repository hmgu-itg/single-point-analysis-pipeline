# Run
1. Modify the `config.yaml` file 
2. Variant QC
3. Single-cohort association
4. Meta-analysis
5. Scan signal regions
6. Plot signals
7. Select independent signals

## 1. Modify the `config.yaml` file
You'll need to provide 5 information.

**input**: Path to the input phenotype files including the `cohort`, `group` and `phenotype` wildcards.  
The format of the phenotype file needs to be as shown below. 
|ID|res|zres|
|---|---|---|
|SampleA|0.5532|0.312|
|SampleB|0.323|0.150|

The third column will be used as the phenotype value.


**cohorts**: For each cohort, the plink binary prefix (`bfile`) and the path to the genetic-relatedness matrix (`grm`) created using GCTA is required.

**group**: An identifier for the phenotype group. This currently has no functional importance other than being used as a prefix for output file names.

**phenotypes**: Either a list of phenotype names or a path to a file containing a list of phenotype names which will be used to fetch the phenotype files based on the `input` configuration value.  

## 2. Variant QC
This step filters the genotype data according to the thresholds specified in the `config.yaml` file.
```bash
snakemake --cores 1 --use-singularity all1
```


## 3. Single-cohort association
Run single-cohort single-point association analysis using GCTA-MLMA
```bash
snakemake --cores 1 --use-singularity all2
```


## 4. Meta-analysis
Run single-point association meta-analysis using METAL
```bash
snakemake --cores 1 --use-singularity create_all_metal
```

## 5. Scan signal regions
Scan the meta-analysis regions to extract significant signal regions.
```bash
snakemake --cores 1 --use-singularity detect_all_peaks
```

## 6. Plot signals
Run PeakPlotter on all signal regions to create regional plots and annotated data.
```bash
snakemake --cores 1 --use-singularity collect_all_peak_csvs
```

## 7. Select independent signals
Run LD-clumping and GCTA-COJO to identify independent signals within all signal regions. 
```bash
snakemake --cores 1 --use-singularity run_all_cojo
```


## Questions
Q. Why do the `freq` and `freq_geno` column values in the `.jma.cojo` file differ?
A. `freq_geno` column is the frequency of the `refA` column allele in the input bfile (you can use `plink --freq` to check).
The `freq` column value is the exact value extracted from the input cojofile, where the cojofile was created from the corresponding metal file.
So the `freq` column value comes from the `Alt_Freq` column value in the metal file, and the `Alt_Freq` column value is the "weighted average of frequency for Alt allele across all studies".
The `freq_geno` and `freq` column values differ because `freq_geno` is just the allele frequency of the variant from the genotype file (plink bfile) that was combined from all cohorts, 
whereas `freq` column is the weighted average of frequency across cohorts (calculated by metal). 

Q. When I try to run a rule, I get an error saying `Text file busy`. What do I do?
A. Delete the script and restore it using `git restore workflow/script/problematic_script.sh`. Your rules should run normally after doing this


## Snakefile order
0. read-config.smk
0. variant-qc.smk
1. single-cohort.smk
2. meta-analysis.smk
3. detect-peaks.smk
4. peakplot.smk
5. cojo.smk
6. query.smk
7. gwas.smk

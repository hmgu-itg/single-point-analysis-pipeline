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



## Questions
Q. Why do the `freq` and `freq_geno` column values in the `.jma.cojo` file differ?
A. `freq_geno` column is the frequency of the `refA` column allele in the input bfile (you can use `plink --freq` to check).
The `freq` column value is the exact value extracted from the input cojofile, where the cojofile was created from the corresponding metal file.
So the `freq` column value comes from the `Alt_Freq` column value in the metal file, and the `Alt_Freq` column value is the "weighted average of frequency for Alt allele across all studies".
The `freq_geno` and `freq` column values differ because `freq_geno` is just the allele frequency of the variant from the genotype file (plink bfile) that was combined from all cohorts, 
whereas `freq` column is the weighted average of frequency across cohorts (calculated by metal). 

Q. When I try to run a rule, I get an error saying `Text file busy`. What do I do?
A. Delete the script and restore it using `git restore workflow/script/problematic_script.sh`. Your rules should run normally after doing this
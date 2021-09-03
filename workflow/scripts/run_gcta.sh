#!/usr/bin/env bash
input_pheno_file=$1
out_prefix=$2
bfile=$3
grm=$4
threads=$5


pheno_file=${out_prefix}.pheno
mlma=${out_prefix}.mlma
mlma_bgz=${out_prefix}.mlma.gz

awk '{OFS="\t"}NR!=1{print $1,$1,$3}' $input_pheno_file > $pheno_file

gcta64 --mlma \
          --bfile $bfile \
          --grm $grm \
          --pheno $pheno_file \
          --out $out_prefix \
          --threads $threads
bgzip -c $mlma > $mlma_bgz
tabix --skip-lines 1 --sequence 1 --begin 3 --end 3 $mlma_bgz

rm $pheno_file $mlma

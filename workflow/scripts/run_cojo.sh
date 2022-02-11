#!/usr/bin/env bash
run_cojo=v0.2.0
bfile=$1
cojofile=$2
metal_file=$3
threshold=$4 # p-value threshold
group=$5
phenotype=$6
chrom=$7
start=$8
end=$9
prefix=${10} # output prefix


echo "=== Running run_cojo.sh ===
run_cojo: $run_cojo
group: $group
phenotype: $phenotype
bfile: $bfile
chrom: $chrom
start: $start
end: $end
threshold: $threshold
cojofile: $cojofile
metal_file: $metal_file
prefix: $prefix
===========================
"
## Subset bfile
echo "[$(date)] Subset bfile"
plink --allow-no-sex --make-bed \
    --bfile $bfile \
    --chr $chrom \
    --from-bp $start \
    --to-bp $end \
    --out $prefix

echo -e "\n\n"
echo "[$(date)] GCTA-COJO"
gcta64 \
  --bfile $prefix \
  --cojo-file <(zcat $cojofile) \
  --out $prefix \
  --cojo-slct \
  --cojo-p $threshold \
  --cojo-wind 10000 \
  --diff-freq 0.2 \
  --cojo-collinear 0.9 || true

echo -e "\n\n"

echo "[$(date)] Done"

#!/usr/bin/env bash
run_cojo=v0.3.0
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
memory=${11}


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
# TODO: Find a way to insert the $threshold variable
# inside the below awk command.

# See Plink 1.9's .qassoc file format for more details
qassoc=$prefix.qassoc
echo "[$(date), Step 1/4] Creating $qassoc"
tabix $metal_file ${chrom}:${start}-${end} | \
    awk 'BEGIN{OFS="\t";print "CHR\tSNP\tBP\tBETA\tSE\tR2\tT\tP"} \
      {
        split($12, a, "[eE]") 
        if($12<1e-6 || (length(a)>1 && a[2] < -6)){
          print $1, "chr"$1":"$3, $3, $10, $11, 1, 1, $12}
      }' > $qassoc


echo
echo "Extracted $(( $(wc -l < $prefix.qassoc) - 1 )) SNP(s)."
echo

## Subset bfile
echo "[$(date), Step 2/4] Subset bfile"
plink --allow-no-sex --make-bed \
    --bfile $bfile \
    --chr $chrom \
    --from-bp $start \
    --to-bp $end \
    --memory $memory \
    --out $prefix

echo -e "\n\n"
echo "[$(date), Step 3/4] Clumping"
## CLUMP from merged file 
plink \
  --bfile $prefix \
  --clump $qassoc \
  --clump-kb 1000 \
  --clump-r2 0.05 \
  --memory $memory \
  --out $prefix

clumpedlist=$prefix.clumped.list
awk -v threshold=$threshold '{if($5<threshold && $1~/^[0-9]+$/){print $3}}' $prefix.clumped > $clumpedlist



echo -e "\n\n"
echo "[$(date), Step 4/4] GCTA-COJO"
gcta64 \
  --bfile $prefix \
  --cojo-file <(zcat $cojofile) \
  --extract $clumpedlist \
  --out $prefix \
  --cojo-slct \
  --cojo-p $threshold \
  --cojo-wind 10000 \
  --diff-freq 0.2 \
  --cojo-collinear 0.9 || true

echo -e "\n\n"

echo "[$(date)] Done"

rm -f \
  $prefix.nosex \
  $prefix.clumped.list
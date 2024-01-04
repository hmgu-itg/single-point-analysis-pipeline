#!/usr/bin/env bash
run_cojo=v0.4.0
bfile=$1
cojofile=$2
assoc_file=$3
threshold=$4 # p-value threshold
chrom=$5
start=$6
end=$7
prefix=${8} # output prefix
memory=${9}


echo "=== Running run_cojo.sh ===
run_cojo: $run_cojo
bfile: $bfile
chrom: $chrom
start: $start
end: $end
threshold: $threshold
cojofile: $cojofile
assoc_file: $assoc_file
prefix: $prefix
===========================
"
# TODO: Create unit test to test whether the threshold values are used correctly
# Determine if the user-defined threshold is less than 1e-6
effective_threshold=$(echo "$threshold" | awk '{if ($1 < 1e-6) print 1e-6; else print $1}')

# See Plink 1.9's .qassoc file format for more details
qassoc=$prefix.qassoc
echo "[$(date), Step 1/4] Creating $qassoc"
tabix $assoc_file ${chrom}:${start}-${end} | \
awk -v threshold_val="$effective_threshold" 'BEGIN{OFS="\t";print "CHR\tSNP\tBP\tBETA\tSE\tR2\tT\tP"}
{
    # Determine if $9 is in scientific notation and extract its exponent
    split($9, a, "[eE]") 
    num_exp = (length(a) > 1) ? a[2] : log($9)/log(10)
    threshold_exp = log(threshold_val)/log(10)

    # Check if $9 is less than the effective threshold
    if($9<threshold_val || num_exp < threshold_exp){
        print $1, $2, $3, $7, $8, 1, 1, $9}
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
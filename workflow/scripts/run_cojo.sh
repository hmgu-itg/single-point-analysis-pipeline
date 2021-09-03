#!/usr/bin/env bash
run_cojo=v0.0.1
bfile=$1
cojofile=$2
metal_file=$3
group=$4
phenotype=$5
chrom=$6
start=$7
end=$8
prefix=$9 # output prefix


echo "=== Running run_cojo.sh ===
run_cojo: $run_cojo
group: $group
phenotype: $phenotype
bfile: $bfile
chrom: $chrom
start: $start
end: $end
cojofile: $cojofile
metal_file: $metal_file
prefix: $prefix
===========================
"

qassoc=$prefix.qassoc
echo "Creating $qassoc"
tabix $metal_file ${chrom}:${start}-${end} | \
    awk 'BEGIN{OFS="\t";print "CHR\tSNP\tBP\tBETA\tSE\tR2\tT\tP"} \
      {if($12<1e-6){ \
          print $1, "chr"$1":"$3, $3, $10, $11, 1, 1, $12} \
      }' > $qassoc
      
      
echo
echo Extracted $(( $(wc -l < $prefix.qassoc) - 1 )) SNPs.
echo

## Subset bfile

plink --allow-no-sex --make-bed \
    --bfile $bfile \
    --chr $chrom \
    --from-bp $start \
    --to-bp $end \
    --out $prefix

echo -e "\n\n"

## CLUMP from merged file 
plink \
  --bfile $prefix \
  --clump $qassoc \
  --clump-kb 1000 \
  --clump-r2 0.05 \
  --out $prefix

clumpedlist=$prefix.clumped.list
awk '$5<5e-8 && $1~/^[0-9]+$/{print $3}' $prefix.clumped > $clumpedlist


filtered_cojofile=$prefix.cojofile
zcat $cojofile | sed 's/:/ /' | awk -v chrom="chr$chrom" -v start="$start" -v end="$end" '{if ($1==chrom && start<=$2 && $2<=end) {print}}' | sed 's/ /:/' > $filtered_cojofile

excludelist=$prefix.excludelist
comm -13 <(cut -f2 $prefix.bim | sort) <(cut -d' ' -f 1 $filtered_cojofile | sort) > $excludelist

if [[ $(wc -l < "$excludelist") -ge 1 ]] ; then
  echo "[WARNING] Variants not present in .bim file but present in cojofile detected"
  echo "[WARNING] Excluding these variants. See list of excluded variants in $excludelist"
  grep -vw -f $excludelist $filtered_cojofile | sponge $filtered_cojofile
fi


echo -e "\n\n"

gcta64 \
  --bfile $prefix \
  --cojo-file $filtered_cojofile \
  --extract $clumpedlist \
  --out $prefix \
  --cojo-slct \
  --cojo-collinear 0.9 || true

bad_freq=$prefix.freq.badsnps
if [[ -f "$bad_freq" ]] ; then
  updated_alleles=$prefix.update_alleles
  flipped=$prefix.merged.flipped
  tail -n+2 $bad_freq | awk -F$'\t' '{print $1,$2,$3,$3,$2}' > $updated_alleles

  echo -e "\n\n"
  plink \
    --bfile $prefix \
    --update-alleles $updated_alleles \
    --out $flipped \
    --make-bed

  echo -e "\n\n"
  gcta64 \
    --bfile $flipped \
    --cojo-file $filtered_cojofile \
    --extract $clumpedlist \
    --out $prefix \
    --cojo-slct \
    --cojo-collinear 0.9
fi

echo -e "\n\n"
echo Done running Cojo

rm \
  $prefix.bed \
  $prefix.bim \
  $prefix.fam \
  $prefix.nosex \
  $prefix.log \
  $prefix.clumped \
  $prefix.clumped.list \
  $prefix.excludelist \
  $prefix.update_alleles \
  $prefix.freq.badsnps \
  $prefix.merged.* \
  $prefix.cojofile
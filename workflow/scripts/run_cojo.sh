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


## CLUMP from merged file 
plink \
  --bfile $prefix \
  --clump $qassoc \
  --clump-kb 1000 \
  --clump-r2 0.05 \
  --out $prefix

clumpedlist=$prefix.clumped.list
awk '$5<5e-8 && $1~/^[0-9]+$/{print $3}' $prefix.clumped > $clumpedlist

##
excludelist=$prefix.excludelist
export RUN_COJO_BIM_VARIABLE="$prefix.bim"
export RUN_COJO_COJOFILE_VARIABLE="$cojofile"
export RUN_COJO_OUTPUTFILE_VARIABLE="$excludelist"
python3 -c "
import os
import pandas as pd
bim = os.environ['RUN_COJO_BIM_VARIABLE']
cojofile = os.environ['RUN_COJO_COJOFILE_VARIABLE']
outputfile = os.environ['RUN_COJO_OUTPUTFILE_VARIABLE']
bim_d = pd.read_csv(bim, sep = '\t', header = None, names = ['chrom', 'SNP', '_', 'pos', 'alt', 'ref'])
cojofile_d = pd.read_csv(cojofile, sep = ' ', header = 0)
cojofile_d = pd.merge(cojofile_d, bim_d, on = 'SNP', how = 'inner')
to_exclude = cojofile_d.loc[cojofile_d['A1']!=cojofile_d['alt'], ['SNP', 'A1', 'A2']]
to_exclude.to_csv(outputfile, sep = ' ', header = False, index = False)
"



filtered_cojofile=$prefix.cojofile # Is this needed??
zcat $cojofile | grep chr$chrom > $filtered_cojofile # Is this needed??

if [[ $(wc -l < "$excludelist") -ge 1 ]] ; then
  echo "[WARNING] Variants not present in .bim file but present in cojofile detected"
  echo "[WARNING] Excluding these variants. See list of excluded variants in $excludelist"
  grep -vw -f $excludelist $filtered_cojofile | sponge $filtered_cojofile
fi


if [ -f "$missnp_file" ] ; then
  echo "[WARNING] excluding multiallelic variants within region in $missnp_file"
  echo "[WARNING] from $filtered_cojofile"
  grep -vw -f $missnp_file $filtered_cojofile | sponge $filtered_cojofile
fi

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
  merged_flipped=$prefix.merged.flipped
  tail -n+2 $bad_freq | awk -F$'\t' '{print $1,$2,$3,$3,$2}' > $updated_alleles
  plink \
    --bfile $merged \
    --update-alleles $updated_alleles \
    --out $merged_flipped \
    --make-bed

  gcta64 \
    --bfile $merged_flipped \
    --cojo-file $filtered_cojofile \
    --extract $clumpedlist \
    --out $prefix \
    --cojo-slct \
    --cojo-collinear 0.9

fi
echo Done running Cojo
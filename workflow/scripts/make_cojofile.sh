#!/usr/bin/env bash
# 
# We need to make sure that we're correctly parsing the effect allele to A1
# and other allele to A2 in the cojofile. METAL outputs Allele1 and Allele2
# as effect and other allele, so we parse Allele1 as A1 and Allele2 as A2. 
# 
# METAL Doc: https://genome.sph.umich.edu/wiki/METAL_Documentation
# GCTA-COJO: https://yanglab.westlake.edu.cn/software/gcta/#COJO


make_cojofile=v0.0.2
group=$1
phenotype=$2
bfile=$3
metal_file=$4
prefix=$5

echo "=== Running make_cojofile.sh ===
run_cojo: $make_cojofile
group: $group
phenotype: $phenotype
bfile: $bfile
metal_file: $metal_file
prefix: $prefix
================================"

samplesize=$prefix.samplesize
plink --bfile $bfile --missing --out $prefix
awk -F '[[:space:]]+' '{if(NR!=1){print $3, $5-$4}}' $prefix.lmiss > $samplesize


tmpfile=$(mktemp /tmp/make_cojofile.XXXXXX) # https://unix.stackexchange.com/a/181938
echo "[INFO] Creating temporary file: $tmpfile"
exec 3>"$tmpfile"

join \
  -j 1 \
  -o 0 1.2 1.3 1.4 1.5 1.6 1.7 2.2 \
  <(zcat $metal_file \
      | awk '{if(NR!=1){print "chr"$1":"$3, $4, $5, $6, $10, $11, $12}}' \
      | sort -k 1b,1 \
   ) \
  <(sort -k 1b,1 $samplesize) >&3
# 'sort -k 1b,1' is the recommended sort command if joining on the first field using 'join'
# See 'man join'

cat \
  <(echo "SNP A1 A2 freq b se p N") \
  <(sed 's/:/ /' $tmpfile \
    | sort -n -k1.4 -k2 \
    | sed 's/ /:/' \
   ) \
  | gzip > $prefix.ma.gz

exec 3>-
# rm $tmpfile
# rm $prefix.{lmiss,samplesize}
rm $prefix.{imiss,nosex}

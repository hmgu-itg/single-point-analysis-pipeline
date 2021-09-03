#!/usr/bin/env bash
make_cojofile=v0.0.1
group=$1
phenotype=$2
bfile=$3
metal_file=$4
output=$5

echo "=== Running make_cojofile.sh ===
run_cojo: $make_cojofile
group: $group
phenotype: $phenotype
bfile: $bfile
metal_file: $metal_file
output: $output
"


plink --bfile $bfile --missing --out $bfile
awk -F '[[:space:]]+' '{if(NR!=1){print $3, $5-$4}}' $bfile.lmiss > $bfile.samplesize

metal_file=CMET.VCAM1.metal.filtered.gz


tmpfile=$(mktemp /tmp/make_cojofile.XXXXXX) # https://unix.stackexchange.com/a/181938
exec 3>"$tmpfile"


join \
  -j 1 \
  -o 0 1.2 1.3 1.4 1.5 1.6 1.7 2.2 \
  <(zcat $metal_file \
      | awk '{if(NR!=1){print "chr"$1":"$3, $5, $4, $6, $10, $11, $12}}' \
      | sort \
   ) \
  <(sort $bfile.samplesize) >&3

cat \
  <(echo "SNP A1 A2 freq b se p N") \
  <(sed 's/:/ /' $tmpfile \
    | sort -n -k1.4 -k2 \
    | sed 's/ /:/' \
    | awk '$4>0 && $NF>500' \
   ) \
  | gzip > $output

exec 3>-
rm $tmpfile
rm $bfile.samplesize $bfile.imiss $bfile.lmiss


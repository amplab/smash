#!/usr/bin/env bash
# Evaluate a VCF on a subset of the genome.

subset=$1 # subset of genome in BED format
true=$2 # true VCF
pred=$3 # predicted VCF
basepath=$4 # path to input data
outpath=$5 # path to output directory 
include_partial=$6 # specifies overlap required for intersection 
shift 6

true_subset=$outpath/${true##*/}.subset
pred_subset=$outpath/${pred##*/}.subset

# -f flag specifies overlap required for intersection
# value of 1.0 requires the entire variant to be in both VCFs
#    to include it in the intersection

function restrict {
    bedtools intersect -a $1 -b $2 -u -f $3
}
restrict $true $subset $include_partial > $true_subset
restrict $pred $subset $include_partial > $pred_subset
subsets="$true_subset $pred_subset"
$basepath/smash/src/python/bench.py $subsets $@
rm $subsets

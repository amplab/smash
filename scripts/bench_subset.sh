#!/usr/bin/env bash

#Copyright (c) 2013, Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice,
#this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice,
#this list of conditions and the following disclaimer in the documentation
#and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

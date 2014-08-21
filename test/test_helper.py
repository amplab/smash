#!/bin/python

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

import StringIO
import vcf
import sys

sys.path.insert(0,'..')
from smashbenchmarking.vcf_eval.variants import Variants
from smashbenchmarking.vcf_eval.chrom_variants import ChromVariants
from smashbenchmarking.parsers.genome import Genome
from smashbenchmarking.normalize_vcf import normalize
from smashbenchmarking.bench import get_contig_lookup

MAX_INDEL_LEN = 50

def str_to_VcfReader(vcf_str):
    str_io = StringIO.StringIO(vcf_str)
    return vcf.Reader(str_io)

def vcf_to_Variants(vcf_str):
    str_io = StringIO.StringIO(vcf_str)
    str_vcf = vcf.Reader(str_io)
    return Variants(str_vcf,MAX_INDEL_LEN)

def vcf_to_ChromVariants(vcf_str,chrom):
    str_io = StringIO.StringIO(vcf_str)
    str_vcf = vcf.Reader(str_io)
    str_vars = Variants(str_vcf,MAX_INDEL_LEN)
    return str_vars.on_chrom(chrom)

def normalize_vcf_to_ChromVariants(vcf_str,chrom):
    str_io = StringIO.StringIO(vcf_str)
    str_vcf = vcf.Reader(str_io)
    norm_iter = normalize(get_reference(),str_vcf)
    str_vars = Variants(norm_iter,MAX_INDEL_LEN)
    return str_vars.on_chrom(chrom)

def get_empty_ChromVariants(contig):
    return ChromVariants(contig,MAX_INDEL_LEN)

def get_reference():
    return Genome('ref.fasta',lambda t: t.split()[0])

def get_reference_index():
    return get_contig_lookup('ref.fasta.fai')
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

import sys
import vcf
import unittest
import StringIO

sys.path.insert(0,'..')
from smashbenchmarking.vcf_eval.rectify_seq import *
from smashbenchmarking.vcf_eval.rectify_seq import _get_chopped_variant,_enlarge_bounds,_get_seq
from smashbenchmarking.vcf_eval.variants import Variants
from smashbenchmarking.parsers.genome import Genome

MAX_INDEL_LEN = 50

def vcf_to_ChromVariants(vcf_str,chrom):
    str_io = StringIO.StringIO(vcf_str)
    str_vcf = vcf.Reader(str_io)
    str_vars = Variants(str_vcf,MAX_INDEL_LEN)
    return str_vars.on_chrom(chrom)

def get_reference():
    return Genome('ref.fasta',lambda t: t.split()[0])

class SequenceRescuerTestCase(unittest.TestCase):
    def testWindowTooBig(self):
        longsv1 = 'ATTGTTCATGA'*300
        longsv2 = 'GCCTAGGGTCA'*300
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   7001   .       """ + longsv1 + """     A       20      PASS    .       GT      1/1\n
chr1   10100   .       """ + longsv2 + """     G       20      PASS    .       GT      0/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   10049   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        fn_vars = vcf_to_ChromVariants(fn_str,'chr1')
        fp_vars = vcf_to_ChromVariants(fp_str,'chr1')
        rescuer = SequenceRescuer('chr1',10049,fn_vars,fp_vars,get_reference(),50)
        self.assertFalse(rescuer.rescued)

    def testEmptyWindow(self):
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   8000   .       G     C       20      PASS    .       GT      1/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   10049   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        fn_vars = vcf_to_ChromVariants(fn_str,'chr1')
        fp_vars = vcf_to_ChromVariants(fp_str,'chr1')
        rescuer = SequenceRescuer('chr1',10049,fn_vars,fp_vars,get_reference(),50)
        self.assertFalse(rescuer.rescued)

    def testTooManyPaths(self):
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   10049   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
chr1   10053   .       TGCGT        T       20      PASS    .       GT      0/1\n
chr1   10055   .       GCTAA        G       20      PASS    .       GT      0/1\n
chr1   10057   .       TA           T       20      PASS    .       GT      1/1\n
chr1   10058   .       GC           G       20      PASS    .       GT      1/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   10025   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
chr1   10028   .       TGCGT        T       20      PASS    .       GT      0/1\n
chr1   10029   .       GCTAA        G       20      PASS    .       GT      0/1\n
chr1   10032   .       TA           T       20      PASS    .       GT      1/1\n
"""
        fn_vars = vcf_to_ChromVariants(fn_str,'chr1')
        fp_vars = vcf_to_ChromVariants(fp_str,'chr1')
        rescuer = SequenceRescuer('chr1',10000,fn_vars,fp_vars,get_reference(),50)
        self.assertFalse(rescuer.rescued)

    def testOnlySnps(self):
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       A     C       20      PASS    .       GT      1/1\n
chr1   7   .       C        T       20      PASS    .       GT      0/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   4   .       A     C       20      PASS    .       GT      1/1\n
"""
        fn_vars = vcf_to_ChromVariants(fn_str,'chr1')
        fp_vars = vcf_to_ChromVariants(fp_str,'chr1')
        rescuer = SequenceRescuer('chr1',3,fn_vars,fp_vars,get_reference(),50)
        self.assertFalse(rescuer.rescued)

    def testFullRescue(self):
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   2   .       TGC     TAT       20      PASS    .       GT      1/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   3   .       G     A       20      PASS    .       GT      1/1\n
chr2   4   .       C     T       20      PASS    .       GT      1/1\n
"""
        fn_vars = vcf_to_ChromVariants(fn_str,'chr2')
        fp_vars = vcf_to_ChromVariants(fp_str,'chr2')
        rescuer = SequenceRescuer('chr2',2,fn_vars,fp_vars,get_reference(),50)
        self.assertTrue(rescuer.rescued)
        self.assertEqual(rescuer.windowsRescued,(0,0))

class RectifySeqHelperMethodsTestCase(unittest.TestCase):
    #find the variant that overlaps the specified location and extends farthest either left or right
    #third arg: False = farthest left, True = farthest right
    def testGetChoppedVariant(self):
        #base case
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88012   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        chopped = _get_chopped_variant(variants,88015,False)
        self.assertEqual(chopped.pos,88012)
        chopped = _get_chopped_variant(variants,88015,True)
        self.assertEqual(chopped.pos,88012)

        #no variants within range
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   87962   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        chopped = _get_chopped_variant(variants,88015,False)
        self.assertFalse(chopped)
        chopped = _get_chopped_variant(variants,88015,True)
        self.assertFalse(chopped)

        #ignore non overlapping snp
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88008   .       T        G       20        PASS     .       GT      1/1\n
chr19   88012   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        chopped = _get_chopped_variant(variants,88015,False)
        self.assertEqual(chopped.pos,88012)
        chopped = _get_chopped_variant(variants,88015,True)
        self.assertEqual(chopped.pos,88012)
        #ignore overlapping snp since indel is farther left (and also right)
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88012   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
chr19   88014   .       T            G       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        chopped = _get_chopped_variant(variants,88015,False)
        self.assertEqual(chopped.pos,88012)
        chopped = _get_chopped_variant(variants,88015,True)
        self.assertEqual(chopped.pos,88012)

        #find longest of overlapping indels
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88008   .       ATTGCTTAACG       A       20       PASS   .       GT      0/1\n
chr19   88012   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        chopped = _get_chopped_variant(variants,88015,False)
        self.assertEqual(chopped.pos,88008)
        chopped = _get_chopped_variant(variants,88015,True)
        self.assertEqual(chopped.pos,88012)
        #test right direction

    def testEnlargeBounds(self):
        #no variant overlaps low, so no change
        #variant exactly abuts high, so get back high + 1
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88012   .       CTTAAGCT     C       20      PASS    .       GT      1/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr19')
        (low,high) = _enlarge_bounds(variants,88000,88020)
        self.assertEqual(low,88000)
        self.assertEqual(high,88021)

    def testGetSeq(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr3   2   .       TCGA     T       20      PASS    .       GT      1/1\n
chr3   9   .       A        AAAA    20      PASS    .       GT      0/1\n
"""
        variants = vcf_to_ChromVariants(pred_str,'chr3')
        window_tup = (1,13,'chr3')
        sequence = _get_seq(window_tup,variants.getAllVariants(),get_reference(),False)
        self.assertEqual(sequence[0],'ATTCGAAAATCG')
        self.assertEqual(sequence[1],'')
        sequence = _get_seq(window_tup,variants.getAllVariants(),get_reference(),True)
        self.assertEqual(sequence[0],'ATTCGATCG')
        self.assertEqual(sequence[1],'ATCGATCGAAAATCG')
if __name__ == "__main__":
    unittest.main()

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
from smashbenchmarking.vcf_eval.eval_helper import *
from smashbenchmarking.vcf_eval.eval_helper import _genotype_concordance_dict
from smashbenchmarking.vcf_eval.chrom_variants import Variant,VARIANT_TYPE,GENOTYPE_TYPE
from smashbenchmarking.vcf_eval.variants import Variants
from smashbenchmarking.parsers.genome import Genome

EPS_BP = 10
EPS_LEN = 10

# these helpers want to live in some common test class
MAX_INDEL_LEN = 50

def vcf_to_ChromVariants(vcf_str,chrom):
    str_io = StringIO.StringIO(vcf_str)
    str_vcf = vcf.Reader(str_io)
    str_vars = Variants(str_vcf,MAX_INDEL_LEN)
    return str_vars.on_chrom(chrom)

def get_reference():
    return Genome('ref.fasta',lambda t: t.split()[0])


class EvalHelperTestCase(unittest.TestCase):
    def testIndelOrSvMatch(self):
        # variant types don't match
        pred_var = Variant(10,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        true_var = Variant(10,'ACGT',['ATGC'],VARIANT_TYPE.INDEL_OTH,GENOTYPE_TYPE.HOM_VAR)
        self.assertFalse(indel_or_sv_match(pred_var,true_var,EPS_BP,EPS_LEN))
        # pred var too far away from true var
        pred_var = Variant(40,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        true_var = Variant(20,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        self.assertFalse(indel_or_sv_match(pred_var,true_var,EPS_BP,EPS_LEN))
        # pred pos off from true pos but within tolerance
        pred_var = Variant(18,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        true_var = Variant(20,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        self.assertTrue(indel_or_sv_match(pred_var,true_var,EPS_BP,EPS_LEN))
        # pred gains too different from true gains
        pred_var = Variant(20,'A',['AAAAAAAAAAAAAAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        true_var = Variant(20,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        self.assertFalse(indel_or_sv_match(pred_var,true_var,EPS_BP,EPS_LEN))
        # pred gains off from true gains but within tolerance
        pred_var = Variant(20,'A',['AAAAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        true_var = Variant(20,'A',['AAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HOM_VAR)
        self.assertTrue(indel_or_sv_match(pred_var,true_var,EPS_BP,EPS_LEN))

    def testRescueMission(self):
        # false negative variant at location is SV; don't rescue
        fn_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   8000   .       G     GATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCT       20      PASS    .       GT      1/1\n
"""
        fp_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   8000   .       G     GC       20      PASS    .       GT      1/1\n
"""
        true_vars = vcf_to_ChromVariants(fn_str,'chr1')
        pred_vars = vcf_to_ChromVariants(fp_str,'chr1')
        num_new_tp,num_removed_fn = rescue_mission(true_vars,pred_vars,8000,get_reference(),100)
        self.assertFalse(any(map(lambda x: x > 0, num_new_tp.itervalues())))
        self.assertFalse(any(map(lambda x: x > 0, num_removed_fn.itervalues())))
        # variant couldn't be rescued; no change to counts or ChromVariants
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
        num_new_tp,num_removed_fn = rescue_mission(fn_vars,fp_vars,2,get_reference(),100)
        self.assertFalse(any(map(lambda x: x > 0, num_new_tp.itervalues())))
        self.assertFalse(any(map(lambda x: x > 0, num_removed_fn.itervalues())))
        self.assertEqual(len(fn_vars.all_locations),2)
        self.assertEqual(len(fp_vars.all_locations),1)
        # variant is rescued; counts change; variants are removed from fn/fp
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
        num_new_tp,num_removed_fn = rescue_mission(fn_vars,fp_vars,2,get_reference(),100)
        self.assertEqual(num_new_tp[VARIANT_TYPE.INDEL_OTH],1)
        self.assertEqual(num_removed_fn[VARIANT_TYPE.SNP],2)
        self.assertEqual(len(fn_vars.all_locations),0)
        self.assertEqual(len(fp_vars.all_locations),0)
    def testChromEvaluateVariants(self):
        # TODO refactor this method
        # TODO test known false positive functionality
        # TODO test genotype concordance
        pass

class ChromVariantStatsTestCase(unittest.TestCase):
    def testInit(self):
        # test counts of false positive, false negative, true positive
        true_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   3   .       G     A       20      PASS    .       GT      1/1\n
chr2   5   .       C     T       20      PASS    .       GT      1/1\n
"""
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   3   .       G     A       20      PASS    .       GT      1/1\n
chr2   7   .       G     C       20      PASS    .       GT      1/1\n
"""
        true_vars = vcf_to_ChromVariants(true_str,'chr2')
        pred_vars = vcf_to_ChromVariants(pred_str,'chr2')
        gtdict = _genotype_concordance_dict()
        gtdict[VARIANT_TYPE.SNP][GENOTYPE_TYPE.HOM_VAR][GENOTYPE_TYPE.HOM_VAR] += 1
        cvs = ChromVariantStats(true_vars,pred_vars,[3],[7],[5],gtdict)
        self.assertEqual(cvs.num_true[VARIANT_TYPE.SNP],2)
        self.assertEqual(cvs.num_pred[VARIANT_TYPE.SNP],2)
        self.assertEqual(len(cvs.false_positives.all_locations),1)
        self.assertEqual(len(cvs.false_negatives.all_locations),1)
        self.assertEqual(cvs.num_tp[VARIANT_TYPE.SNP],1)
        self.assertEqual(cvs.num_fn[VARIANT_TYPE.SNP],1)
        self.assertEqual(cvs.num_fp[VARIANT_TYPE.SNP],1)

    def testRectify(self):
        # rectify CVS with a rescue-able indel
        true_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   2   .       TGC     TAT       20      PASS    .       GT      1/1\n
"""
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2   3   .       G     A       20      PASS    .       GT      1/1\n
chr2   4   .       C     T       20      PASS    .       GT      1/1\n
"""
        true_vars = vcf_to_ChromVariants(true_str,'chr2')
        pred_vars = vcf_to_ChromVariants(pred_str,'chr2')
        gtdict = _genotype_concordance_dict() # leave empty, we aren't testing this yet
        cvs = ChromVariantStats(true_vars,pred_vars,[],[3,4],[2],gtdict)
        # before rectify, no true positives
        self.assertTrue(all(map(lambda x: x == 0,cvs.num_tp.itervalues())))
        # one false negative indel
        self.assertEqual(cvs.num_fn[VARIANT_TYPE.INDEL_OTH],1)
        # two false positives SNPs
        self.assertEqual(cvs.num_fp[VARIANT_TYPE.SNP],2)
        cvs.rectify(get_reference(),100)
        # after rectify, one true positive indel
        self.assertEqual(cvs.num_tp[VARIANT_TYPE.INDEL_OTH],1)
        # no false positives or false negatives
        self.assertTrue(all(map(lambda x: x == 0, cvs.num_fp.itervalues())))
        self.assertTrue(all(map(lambda x: x ==0, cvs.num_fn.itervalues())))





if __name__ == "__main__":
    unittest.main()
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

from test_helper import MAX_INDEL_LEN,vcf_to_Variants,get_reference

sys.path.insert(0,'..')
from smashbenchmarking import Variants,evaluate_variants
from smashbenchmarking import VARIANT_TYPE

sv_eps = 100

class EvalVarsTestCase(unittest.TestCase):
    def truePositive(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],1)
        self.assertEqual(stats['num_pred'],1)
        self.assertEqual(stats['good_predictions'],1)
        self.assertEqual(stats['intersect_bad'],0)
        self.assertEqual(stats['false_negatives'],0)
        self.assertEqual(stats['nrd_total'],1)
        self.assertEqual(stats['nrd_wrong'],0)

    def truePositiveBadGeno(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],1)
        self.assertEqual(stats['num_pred'],1)
        self.assertEqual(stats['good_predictions'],1)
        self.assertEqual(stats['intersect_bad'],0)
        self.assertEqual(stats['false_negatives'],0)
        self.assertEqual(stats['nrd_total'],1)
        self.assertEqual(stats['nrd_wrong'],1)

    def falseNegative(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],1)
        self.assertEqual(stats['num_pred'],0)
        self.assertEqual(stats['good_predictions'],0)
        self.assertEqual(stats['intersect_bad'],0)
        self.assertEqual(stats['false_negatives'],1)
        self.assertEqual(stats['nrd_total'],0)
        self.assertEqual(stats['nrd_wrong'],0)

    def falsePositive(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],0)
        self.assertEqual(stats['num_pred'],1)
        self.assertEqual(stats['good_predictions'],1)
        self.assertEqual(stats['intersect_bad'],1)
        self.assertEqual(stats['false_negatives'],1)
        self.assertEqual(stats['nrd_total'],0)
        self.assertEqual(stats['nrd_wrong'],0)

    def trueNegative(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],0)
        self.assertEqual(stats['num_pred'],0)
        self.assertEqual(stats['good_predictions'],0)
        self.assertEqual(stats['intersect_bad'],0)
        self.assertEqual(stats['false_negatives'],0)
        self.assertEqual(stats['nrd_total'],0)
        self.assertEqual(stats['nrd_wrong'],0)

    def badCallAtTrueSite(self,stat_reporter,var_type):
        stats = stat_reporter(var_type)
        self.assertEqual(stats['num_true'],1)
        self.assertEqual(stats['num_pred'],1)
        self.assertEqual(stats['good_predictions'],0)
        self.assertEqual(stats['intersect_bad'],1)
        self.assertEqual(stats['false_negatives'],1)
        self.assertEqual(stats['nrd_total'],0)
        self.assertEqual(stats['nrd_wrong'],0)

    def setUp(self):
        true_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""
        true_io = StringIO.StringIO(true_str)
        true_vcf = vcf.Reader(true_io)
        self.true_vars = Variants(true_vcf, MAX_INDEL_LEN)

    def tearDown(self):
        pass

    def test_perfect(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.truePositive(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.truePositive(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_false_neg_snp(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.falseNegative(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.truePositive(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_bad_snp(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       G       20      PASS    .       GT      0/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.badCallAtTrueSite(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.truePositive(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_false_neg_sv_ins(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        sv_eps = 100

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.truePositive(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.falseNegative(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_sv_snp_collision(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1\n
chr19   269751  .       A       T       20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        sv_eps = 100

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        snp_stats = stat_reporter(VARIANT_TYPE.SNP)

        self.assertEqual(snp_stats['num_true'],1)
        self.assertEqual(snp_stats['num_pred'],2)
        self.assertEqual(snp_stats['good_predictions'],1)
        self.assertEqual(snp_stats['intersect_bad'],0)
        self.assertEqual(snp_stats['false_negatives'],0)

        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)

        sv_ins_stats = stat_reporter(VARIANT_TYPE.SV_INS)

        self.assertEqual(sv_ins_stats['num_true'],1)
        self.assertEqual(sv_ins_stats['num_pred'],0)
        self.assertEqual(sv_ins_stats['good_predictions'],0)
        self.assertEqual(sv_ins_stats['intersect_bad'],1)
        self.assertEqual(sv_ins_stats['false_negatives'],1)

        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)



    def test_snp_sv_collision(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      0/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        sv_eps = 100

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        snp_stats = stat_reporter(VARIANT_TYPE.SNP)

        self.assertEqual(snp_stats['num_true'],1)
        self.assertEqual(snp_stats['num_pred'],0)
        self.assertEqual(snp_stats['good_predictions'],0)
        self.assertEqual(snp_stats['intersect_bad'],1)
        self.assertEqual(snp_stats['false_negatives'],1)

        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)

        sv_ins_stats = stat_reporter(VARIANT_TYPE.SV_INS)

        self.assertEqual(sv_ins_stats['num_true'],1)
        self.assertEqual(sv_ins_stats['num_pred'],2)
        self.assertEqual(sv_ins_stats['good_predictions'],1)
        self.assertEqual(sv_ins_stats['intersect_bad'],0)
        self.assertEqual(sv_ins_stats['false_negatives'],0)

        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_perfect_bad_genotypes(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      1/1\n
chr19   89272   .       C       T       20      PASS    .       GT      1/1\n
chr19   269751  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      0/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.truePositiveBadGeno(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositiveBadGeno(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.truePositiveBadGeno(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_approx_sv(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1\n
chr19   269771  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.truePositive(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)
        self.truePositive(stat_reporter,VARIANT_TYPE.SV_INS)
        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_sv_out_of_range(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   88013   .       CTT     C       20      PASS    .       GT      0/1\n
chr19   89272   .       C       T       20      PASS    .       GT      0/1\n
chr19   269852  .       A       AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG 20      PASS    .       GT      1/1
"""

        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        stat_reporter, errors = evaluate_variants(self.true_vars, pred_vars, sv_eps, sv_eps, None, None, None)

        self.truePositive(stat_reporter,VARIANT_TYPE.SNP)
        self.trueNegative(stat_reporter,VARIANT_TYPE.INDEL_INS)
        self.truePositive(stat_reporter,VARIANT_TYPE.INDEL_DEL)

        sv_ins_stats = stat_reporter(VARIANT_TYPE.SV_INS)

        self.assertEqual(sv_ins_stats['num_true'],1)
        self.assertEqual(sv_ins_stats['num_pred'],1)
        self.assertEqual(sv_ins_stats['good_predictions'],0)
        self.assertEqual(sv_ins_stats['intersect_bad'],0)
        self.assertEqual(sv_ins_stats['false_negatives'],1)
        self.assertEqual(sv_ins_stats['nrd_total'],0)
        self.assertEqual(sv_ins_stats['nrd_wrong'],0)

        self.trueNegative(stat_reporter,VARIANT_TYPE.SV_DEL)

    def test_known_false_positives(self):
        true_vcf = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1    1       .       T       A       20      PASS    .       GT       0/1\n
chr1    8       .       A       C       20      PASS    .       GT       1/1\n
"""
        pred_vcf = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1    3       .       G       C       20      PASS     .      GT      1/1\n
chr1    5       .       C       G       20      PASS     .      GT      0/1\n
chr1    8       .       A       C       20      PASS     .      GT      1/1\n
"""
        known_fp_vcf = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1    3       .       G       .       20      PASS    .       GT      0/0\n
chr1    5       .       C       G       20      PASS    .       GT      0/0\n
chr1    9       .       T       .       20      PASS    .       GT      0/0\n
"""
        known_fp_io = StringIO.StringIO(known_fp_vcf)
        known_fp_vars = Variants(vcf.Reader(known_fp_io),MAX_INDEL_LEN,knownFP=True)

        stat_reporter, vcf_output = evaluate_variants(vcf_to_Variants(true_vcf),vcf_to_Variants(pred_vcf),sv_eps,sv_eps, \
            get_reference(),50,known_fp_vars)

        snp_stats = stat_reporter(VARIANT_TYPE.SNP)

        self.assertEqual(snp_stats['num_true'],2)
        self.assertEqual(snp_stats['num_pred'],3)
        self.assertEqual(snp_stats['good_predictions'],1)
        self.assertEqual(snp_stats['false_positives'],2) # predicted vars not in ground truth
        self.assertEqual(snp_stats['false_negatives'],1)
        self.assertEqual(snp_stats['known_fp_calls'],2)
        self.assertEqual(snp_stats['known_fp'],2)

if __name__ == '__main__':
    unittest.main()
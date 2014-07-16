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

from test_helper import MAX_INDEL_LEN

sys.path.insert(0,'..')
from smashbenchmarking.vcf_eval.chrom_variants import Variant,VARIANT_TYPE,GENOTYPE_TYPE
from smashbenchmarking.vcf_eval.chrom_variants import *
from smashbenchmarking.vcf_eval.chrom_variants import _getOverlaps, _getRestOfPath
from smashbenchmarking.vcf_eval.variants import Variants

#this class holds info from each VCF record
class VariantTestCase(unittest.TestCase):

    def testSnpVariant(self):
        testVar = Variant(10,'A',['C'],VARIANT_TYPE.SNP,GENOTYPE_TYPE.HET)
        self.assertEqual(testVar.gains,[0])
        self.assertEqual(testVar.losses,[0])

    def testInsertionVariant(self):
        testVar = Variant(10,'A',['AAAA'],VARIANT_TYPE.INDEL_INS,GENOTYPE_TYPE.HET)
        self.assertEqual(testVar.gains,[3])
        self.assertEqual(testVar.losses,[0])
        self.assertFalse(testVar.overlaps_allele(9))
        self.assertTrue(testVar.overlaps_allele(10))
        self.assertFalse(testVar.overlaps_allele(11))

    def testDeletionVariant(self):
        testVar = Variant(10,'AAAA',['A'],VARIANT_TYPE.INDEL_DEL,GENOTYPE_TYPE.HET)
        self.assertEqual(testVar.gains,[-3])
        self.assertEqual(testVar.losses,[3])
        self.assertFalse(testVar.overlaps_allele(9))
        self.assertTrue(testVar.overlaps_allele(10))
        self.assertTrue(testVar.overlaps_allele(13))
        self.assertFalse(testVar.overlaps_allele(14))

    def testOverlapsVariant(self):
        snpVar = Variant(10,'A',['C'],VARIANT_TYPE.SNP,GENOTYPE_TYPE.HET)
        indelVar = Variant(7,'AAAAAAAA',['A'],VARIANT_TYPE.INDEL_DEL,GENOTYPE_TYPE.HET)
        self.assertTrue(snpVar.strictly_overlaps_var(indelVar))
        self.assertTrue(indelVar.strictly_overlaps_var(snpVar))
        snpVar2 = Variant(2,'C',['T'],VARIANT_TYPE.SNP,GENOTYPE_TYPE.HET)
        self.assertFalse(snpVar2.strictly_overlaps_var(indelVar))
        self.assertFalse(indelVar.strictly_overlaps_var(snpVar2))

#test ChromVariants class
class ChromVariantsTestCase(unittest.TestCase):

    def testAddRecord(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001
chr19   2       .       A       T       20      PASS    .       GT      0/1\n
chr19   5       .       AT      A       20      PASS    .       GT      1/1\n
chr19   10      .       C       CG      20      PASS    .       GT      0/1\n
chr19   20      .       ATGC    ACGT    20      PASS    .       GT      0/1\n
chr19   30  .       AAAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG       A 20      PASS    .       GT      0/1
"""
        newChromVar = ChromVariants('chr19',MAX_INDEL_LEN)
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        newChromVar.add_record(pred_vcf.next())
        self.assertTrue(2 in newChromVar.all_locations)
        self.assertEqual(newChromVar.snp_pos_dict[2].var_type,VARIANT_TYPE.SNP)
        self.assertFalse(newChromVar.indel_pos_dict)
        self.assertFalse(newChromVar.sv_pos_dict)
        newChromVar.add_record(pred_vcf.next())
        self.assertTrue(5 in newChromVar.all_locations)
        self.assertEqual(newChromVar.indel_pos_dict[5].var_type, VARIANT_TYPE.INDEL_DEL)
        newChromVar.add_record(pred_vcf.next())
        self.assertTrue(10 in newChromVar.all_locations)
        self.assertEqual(newChromVar.indel_pos_dict[10].var_type, VARIANT_TYPE.INDEL_INS)
        newChromVar.add_record(pred_vcf.next())
        self.assertTrue(20 in newChromVar.all_locations)
        self.assertEqual(newChromVar.indel_pos_dict[20].var_type, VARIANT_TYPE.INDEL_INV)
        newChromVar.add_record(pred_vcf.next())
        self.assertTrue(30 in newChromVar.all_locations)
        self.assertEqual(newChromVar.sv_pos_dict[30].var_type, VARIANT_TYPE.SV_DEL)
        #all indels/sv (del,ins,oth) live in one bucket
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.INDEL_DEL)),len(newChromVar._var_dict(VARIANT_TYPE.INDEL_INS)))
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.SV_DEL)),len(newChromVar._var_dict(VARIANT_TYPE.SV_INS)))

    def testAddRecordNoSample(self):
        vcf_str = """##fileformat=VCFv4.0\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    \n
chr19   2       .       A       T       20      PASS    .       \n
"""
        newChromVar = ChromVariants('chr19',MAX_INDEL_LEN)
        test_vcf = vcf.Reader(StringIO.StringIO(vcf_str))
        newChromVar.add_record(test_vcf.next())
        self.assertEqual(newChromVar.all_locations,[])

    def testRemoveRecord(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001
chr19   1       .       C       G       20      PASS    .       GT      1/1\n
chr19   2       .       A       T       20      PASS    .       GT      0/1\n
chr19   5       .       AT      A       20      PASS    .       GT      1/1\n
chr19   10      .       C       CG      20      PASS    .       GT      0/1\n
chr19   20      .       ATGC    ACGT    20      PASS    .       GT      0/1\n
chr19   30  .       AAAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG       A 20      PASS    .       GT      0/1
"""
        newChromVar = ChromVariants('chr19',MAX_INDEL_LEN)
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        for r in pred_vcf:
            newChromVar.add_record(r)
        self.assertEqual(len(newChromVar.all_variants),6)
        self.assertEqual(len(newChromVar.all_locations),6)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.SNP]),2)
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.SNP)),2)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.INDEL_DEL]),1)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.INDEL_INS]),1)
        # note that _var_dict holds all types of indels in one bucket, etc
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.INDEL_INS)),3)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.SV_DEL]),1)
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.SV_DEL)),1)

        newChromVar._remove_variant(10)
        self.assertEqual(len(newChromVar.all_variants),5)
        self.assertEqual(len(newChromVar.all_locations),5)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.SNP]),2)
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.SNP)),2)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.INDEL_DEL]),1)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.INDEL_INS]),0)
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.INDEL_INS)),2)
        self.assertEqual(len(newChromVar._var_locations[VARIANT_TYPE.SV_DEL]),1)
        self.assertEqual(len(newChromVar._var_dict(VARIANT_TYPE.SV_DEL)),1)

    def testKnownFalsePositives(self):
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1    7       .       G       .       20      PASS    .       GT       0/0\n
"""
        newChromVar = ChromVariants('chr1',MAX_INDEL_LEN,knownFP=True)
        vcf_io = StringIO.StringIO(vcf_str)
        vcfr = vcf.Reader(vcf_io)
        for r in vcfr:
            newChromVar.add_record(r)
        self.assertEqual(newChromVar.all_locations,[7])
        var = newChromVar.all_variants[7]
        self.assertEqual(var.ref,'G')
        self.assertEqual(var.alt[0], 'NONE')

#test rando helper methods
class ChromVariantHelperMethodsTestCase(unittest.TestCase):
    def testExtractRange(self):
        new_range = extract_range([8,9,10,11,20,80],10,20)
        self.assertEqual(new_range[0],10)
        self.assertEqual(new_range[-1],11)

    def testExtractRangeAndFilter(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   2    .       T     G       20      PASS    .       GT      0/1\n
chr19   10   .       A     C       20      PASS    .       GT      1/1\n
chr19   13   .       A       ACT       20      PASS    .       GT      1/1\n
chr19   15   .       A       T         20      PASS    .       GT      0/1\n
chr19   18  .       AAAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG       A 20      PASS    .       GT      0/1
"""
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)

        variants_in_window = extract_range_and_filter(pred_vars.on_chrom('chr19'),10,20,13)
        self.assertEqual(len(variants_in_window),3)
        #SV is removed
        self.assertFalse(any(map(lambda v: v.var_type.startswith("SV"), variants_in_window)))
        #variant overlapping with variant at location of interest is removed
        self.assertFalse(any(map(lambda v: v.pos == 2,variants_in_window)))

    def testGetOverlaps(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   10   .       ACT     A       20      PASS    .       GT      1/1\n
chr19   13   .       AC      A       20       PASS    .      GT      1/1\n
chr19   14   .       TAGG      TA        20      PASS    .       GT      1/1\n
chr19   15   .       AGG       A         20      PASS    .       GT      0/1\n
chr19   19  .       T       TAAAC 20      PASS    .       GT      0/1
"""
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf,MAX_INDEL_LEN)
        variants_in_window = extract_range_and_filter(pred_vars.on_chrom('chr19'),10,20,10)

        #the three overlapping variants should be in same group
        overlaps = _getOverlaps([], variants_in_window)
        self.assertEqual(len(overlaps),3)
        self.assertEqual(map(lambda o: len(o),overlaps),[1,3,1])

    def testGetRestOfPaths(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   11   .       ACT     A       20      PASS    .       GT      1/1\n
chr19   15   .       ACGATT      AA       20       PASS    .      GT      1/1\n
chr19   16   .       ACG      A        20      PASS    .       GT      1/1\n
chr19   22   .       ATT       A         20      PASS    .       GT      0/1\n
"""
        pred_io = StringIO.StringIO(pred_str)
        pred_vcf = vcf.Reader(pred_io)
        pred_vars = Variants(pred_vcf,MAX_INDEL_LEN)
        viw = extract_range_and_filter(pred_vars.on_chrom('chr19'),10,25,11)

        paths = _getRestOfPath([], _getOverlaps([],viw))
        #all paths take variants at pos 11 and 22; one takes pos 15, one pos 16
        self.assertEqual(len(paths),2)
        self.assertEqual(len(paths[0]),3)
        self.assertEqual(len(paths[1]),3)
        self.assertTrue(all(map(lambda e: any(map(lambda x: x.pos == 11, e)), paths)))
        self.assertTrue(all(map(lambda e: any(map(lambda x: x.pos == 22, e)), paths)))
        self.assertTrue(any(map(lambda x: x.pos == 15, paths[0])))
        self.assertFalse(any(map(lambda x: x.pos == 16, paths[0])))
        self.assertFalse(any(map(lambda x: x.pos == 15, paths[1])))
        self.assertTrue(any(map(lambda x: x.pos == 16, paths[1])))

    def testIsInversion(self):
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   11      .       ACT     TCA       20      PASS    .       GT      1/1\n
chr19   15      .       ACGATT  ATTAGC      20      PASS    .       GT      1/1\n
chr19   16      .       ACG     A       20      PASS    .       GT      1/1\n
"""
        vcfr = vcf_file = vcf.Reader(StringIO.StringIO(pred_str))
        self.assertTrue(is_inversion(vcfr.next(),MAX_INDEL_LEN)) # inversion with no leading base
        self.assertTrue(is_inversion(vcfr.next(),MAX_INDEL_LEN)) # inversion with leading base
        self.assertFalse(is_inversion(vcfr.next(),MAX_INDEL_LEN)) # deletions are not inversions


if __name__ == '__main__':
    unittest.main()

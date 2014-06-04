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

import sys, unittest, vcf
import StringIO

from test_helper import get_reference, vcf_to_ChromVariants, MAX_INDEL_LEN

sys.path.insert(0,'..')
from smashbenchmarking.normalize_vcf import normalize,write,find_redundancy, left_normalize
from smashbenchmarking.parsers.vcfwriter import VCFWriter

class NormalizeTestCase(unittest.TestCase):
    def setUp(self):
        left_slides = []
        self.test_fasta = 'ref.fasta'

    def tearDown(self):
        pass

    def countRecords(self, vcf):
        count = 0
        for r in vcf:
            count += 1
        return count

    def getVcf(self,str):
        vcf_io = StringIO.StringIO(str)
        return vcf.Reader(vcf_io)

    def normalizeStringToWriter(self,vcf_str):
        vcf_io = StringIO.StringIO(vcf_str)
        test_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter('ref.fasta','name',output_io)
        map(lambda r: write(r,output_writer),normalize(get_reference(),test_vcf))
        outputStr = output_io.getvalue()
        outputStr = outputStr.replace('\n','\n\n')
        return vcf.Reader(StringIO.StringIO(outputStr))

    def testFindRedundancy(self):
        alleles = ['acc','gcc','tcc']
        self.assertEqual(find_redundancy(alleles),2)
        alleles = ['atcc','cgtcc']
        self.assertEqual(find_redundancy(alleles),3)
        alleles = ['tcc','atcc']
        self.assertEqual(find_redundancy(alleles),2)
        alleles = ['atg','gta']
        self.assertEqual(find_redundancy(alleles),0)

    def testLeftNormalize(self):
        #left normalize deletion
        norm_pos, norm_ref, norm_alts = left_normalize(get_reference(),'chr1',2,'CGCCG',['CG'])
        self.assertEqual(norm_pos,0)
        self.assertEqual(norm_ref,'AACGC')
        self.assertEqual(norm_alts[0],'AA')

        #left normalize insertion
        norm_pos, norm_ref, norm_alts = left_normalize(get_reference(),'chr4',12,'G',['GGG'])
        self.assertEqual(norm_pos,7)
        self.assertEqual(norm_ref,'C')
        self.assertEqual(norm_alts[0],'CGG')

    def testCleanOnly(self):
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2    6       .       g       cg       20      PASS    .       GT      0/1\n
"""
        norm = normalize(get_reference(),self.getVcf(vcf_str),50,True)
        record = norm.next()
        self.assertEqual(record.POS,6)
        self.assertEqual(record.REF,'G')
        self.assertEqual(record.ALT,['CG'])

    def testNormalize(self):
        #regular records are unchanged
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     A       20      PASS    .       GT      0/1\n
"""
        norm_vcf = normalize(get_reference(),self.getVcf(vcf_str))
        self.assertEqual(self.countRecords(norm_vcf),1)

        #test that hom ref records are removed
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     C       20      PASS    .       GT      0/0\n
chr1   3   .       G     A       20      PASS    .       GT      1/1\n
"""
        norm_vcf = normalize(get_reference(),self.getVcf(vcf_str))
        self.assertEqual(self.countRecords(norm_vcf),1)

        #test that SNP/indels without genotyping are removed
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     A       20      PASS    .       GT      .\n
chr1   3   .       G     C       20      PASS    .       GT      0/0\n
chr1   4   .       G     T       20      PASS    .       GT      0|0\n
chr1   5   .       G     A       20      PASS    .       GT      1/1\n
"""
        norm_vcf = normalize(get_reference(),self.getVcf(vcf_str))
        self.assertEqual(self.countRecords(norm_vcf),1)

        #test that SV without genotyping is retained
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG       20      PASS    .       GT      .\n
"""
        norm_vcf = normalize(get_reference(),self.getVcf(vcf_str))
        self.assertEqual(self.countRecords(norm_vcf),1)

        #test that lower case ref/alt gets upper-cased
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       c     a       20      PASS    .       GT      0/1\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        lowercase_vcf = vcf.Reader(StringIO.StringIO(vcf_str))
        output_vcf = normalize(get_reference(),self.getVcf(vcf_str))
        original_r = lowercase_vcf.next()
        norm_r = output_vcf.next()
        self.assertEqual(original_r.REF,'c')
        self.assertEqual(original_r.ALT[0], 'a')
        self.assertEqual(norm_r.REF,'C')
        self.assertEqual(norm_r.ALT[0],'A')

        # test normalizing an insertion
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   9   .       a     ga       20      PASS    .       GT      0/1\n
"""
        record = normalize(get_reference(),self.getVcf(vcf_str)).next()
        self.assertEqual(record.POS,6)
        self.assertEqual(record.REF,'C')
        self.assertEqual(record.ALT,['CG'])

        # test normalizing a deletion
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   5   .       cc     c       20      PASS    .       GT      0/1\n
"""
        record = normalize(get_reference(),self.getVcf(vcf_str)).next()
        self.assertEqual(record.POS,4)
        self.assertEqual(record.REF,'GC')
        self.assertEqual(record.ALT,['G'])

    def testGenotypes(self):
        # keep genotype info for a compound heterozygous call
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       A     C,T       20      PASS    .       GT      1/2\n
"""
        vcf = self.getVcf(vcf_str)
        record = normalize(get_reference(),vcf).next()
        self.assertEqual(record.samples[0].gt_nums, "1/2")

    def testMultipleAltAlleles(self):
        # multiple alleles aren't normalized if the two alt alleles would be normalized differently
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2    6       .       G       CG       20      PASS    .       GT      0/1\n
"""
        record = normalize(get_reference(),self.getVcf(vcf_str)).next()
        self.assertEqual(record.POS,3)
        self.assertEqual(record.REF,'G')
        self.assertEqual(record.ALT[0], 'GC')
        vcf_str2 = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2    6       .       G       CG,C       20      PASS    .       GT      0/1\n
"""
        record = normalize(get_reference(),self.getVcf(vcf_str2)).next()
        self.assertEqual(record.POS,6)
        self.assertEqual(record.REF,'G')
        self.assertEqual(record.ALT[0],'CG')

    def testNormalizerWriter(self):
        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       a     c       20      PASS    .       GT      0/1\n
chr1   4   .      A      G       20      PASS     .      GT      1/1\n
"""
        output_vcf = self.normalizeStringToWriter(vcf_str)
        r1 = output_vcf.next()
        self.assertEqual(r1.POS,2)

    def testCollidingVariants(self):
        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   5   .      A     TGC       20      PASS    .       GT      1/1\n
chr1   5   .      A      GGG       20      PASS     .      GT      1/1\n
"""
        norm_iter = normalize(get_reference(),self.getVcf(vcf_str))
        count = self.countRecords(norm_iter)
        self.assertEqual(count,1)

    def testNormalizedToCollision(self):
        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr2    4       .       C       T       20      PASS    .       GT      0/1\n
chr2    5       .       C       CGC     20      PASS    .       GT      0/1\n
chr4    2       .       A       AGG     20      PASS    .       GT      0/1\n
chr4    6       .       C       CTC     20      PASS    .       GT      0/1\n
"""
        norm_iter = normalize(get_reference(),self.getVcf(vcf_str))
        r1 = norm_iter.next()
        r2 = norm_iter.next()
        r3 = norm_iter.next()
        r4 = norm_iter.next()
        self.assertEqual(r1.POS,4) # chr2 SNP doesn't change
        self.assertEqual(r2.POS,5) # chr2 insertion gets normed forward 1 base and slid back to original pos
        self.assertEqual(r2.REF,"C")
        self.assertEqual(r2.ALT,["CGC"])
        self.assertEqual(r3.POS,2)
        self.assertEqual(r3.REF,"A")
        self.assertEqual(r3.ALT,["AGG"])
        self.assertEqual(r4.POS,3)
        self.assertEqual(r4.REF,"T")
        self.assertEqual(r4.ALT,["TCT"])

        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr4    2       .       ATC     A     20      PASS    .       GT      0/1\n
chr4    6       .       CTC     C     20      PASS    .       GT      0/1\n
"""
        norm_iter = normalize(get_reference(),self.getVcf(vcf_str))
        r1 = norm_iter.next()
        r2 = norm_iter.next()
        self.assertEqual(r1.POS,2)
        self.assertEqual(r1.REF,"ATC")
        self.assertEqual(r1.ALT,["A"])
        self.assertEqual(r2.POS,5)
        self.assertEqual(r2.REF,"TCT")
        self.assertEqual(r2.ALT,["T"])


    def testNormalizeTwoToCollision(self):
        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr4    4       .       C       CTC     20      PASS    .       GT      0/1\n
chr4    6       .       C       CTC     20      PASS    .       GT      0/1\n
"""
        norm_iter = normalize(get_reference(),self.getVcf(vcf_str))
        r1 = norm_iter.next()
        r2 = norm_iter.next()
        self.assertEqual(r1.POS,2)
        self.assertEqual(r1.REF,"A")
        self.assertEqual(r1.ALT,["ATC"])
        self.assertEqual(r2.POS,3)
        self.assertEqual(r2.REF,"T")
        self.assertEqual(r2.ALT,["TCT"])

    def testNormalizeThreeCollision(self):
        # the OP info flag is fake
        vcf_str = """##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr4    3       .       T       C     20      PASS    .       GT      0/1\n
chr4    3       .       T       TCTTC     20      PASS    OP=1       GT      0/1\n
chr4    3       .       TCTC    T        20      PASS     OP=2        GT     0/1\n
"""
        norm_iter = normalize(get_reference(),self.getVcf(vcf_str))
        r1 = norm_iter.next()
        r2 = norm_iter.next()
        r3 = norm_iter.next()
        self.assertEqual(r1.POS,3)
        self.assertEqual(r2.POS,4)
        self.assertEqual(r2.REF,"C")
        self.assertEqual(r2.ALT,["CTTCC"])
        #self.assertEqual(r3.POS,5)

if __name__ == '__main__':
    unittest.main()

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
from smashbenchmarking.normalize_vcf import find_redundancy, left_normalize, normalize
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

    def outputToVcf(self, outputIO):
        outputStr = outputIO.getvalue()
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
        norm_pos, norm_ref, norm_alts = left_normalize(get_reference(),'chr2',4,'CGGA',['CTTGGA'])
        self.assertEqual(norm_pos,1)
        self.assertEqual(norm_ref,'TGCC')
        self.assertEqual(norm_alts[0],'TGCCTT')

    def testNormalize(self):
        #regular records are unchanged
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     A       20      PASS    .       GT      0/1\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        norm_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,norm_vcf,output_writer)
        self.assertEqual(self.countRecords(self.outputToVcf(output_io)),1)

        #test that hom ref records are removed
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     C       20      PASS    .       GT      0/0\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        homref_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,homref_vcf,output_writer)
        self.assertEqual(self.countRecords(self.outputToVcf(output_io)),0)

        #test that SNP/indels without genotyping are removed
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     A       20      PASS    .       GT      .\n
chr1   3   .       G     C       20      PASS    .       GT      0/0\n
chr1   4   .       G     T       20      PASS    .       GT      0|0\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        homref_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,homref_vcf,output_writer)
        self.assertEqual(self.countRecords(self.outputToVcf(output_io)),0)

        #test that SV without genotyping is retained
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       C     AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG       20      PASS    .       GT      .\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        homref_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,homref_vcf,output_writer)
        self.assertEqual(self.countRecords(self.outputToVcf(output_io)),1)

        #test that lower case ref/alt gets upper-cased
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       c     a       20      PASS    .       GT      0/1\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        lowercase_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,lowercase_vcf,output_writer)
        output_vcf = self.outputToVcf(output_io)
        lowercase_vcf = vcf.Reader(StringIO.StringIO(vcf_str))
        original_r = lowercase_vcf.next()
        norm_r = output_vcf.next()
        self.assertEqual(original_r.REF,'c')
        self.assertEqual(original_r.ALT[0], 'a')
        self.assertEqual(norm_r.REF,'C')
        self.assertEqual(norm_r.ALT[0],'A')

    def testGenotypes(self):
        # normalize a compound heterozygous call
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
##source=TVsim\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1   2   .       A     C       20      PASS    .       GT      1/2\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        test_vcf = vcf.Reader(vcf_io)
        output_io = StringIO.StringIO()
        output_writer = VCFWriter(self.test_fasta,'name',output_io)
        normalize(self.test_fasta,test_vcf,output_writer)
        output_vcf = self.outputToVcf(output_io)
        record = output_vcf.next()
        self.assertEqual(record.samples[0].gt_nums, "1/2")

if __name__ == '__main__':
    unittest.main()

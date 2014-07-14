#!/usr/bin/env python

import os
import sys
import vcf
import unittest
import StringIO

from test_helper import MAX_INDEL_LEN,vcf_to_ChromVariants,get_reference

sys.path.insert(0,'..')
from smashbenchmarking.vcf_eval.variants import *
from smashbenchmarking.vcf_eval.variants import _aggregate
from smashbenchmarking.vcf_eval.chrom_variants import VARIANT_TYPE,GENOTYPE_TYPE
from smashbenchmarking.vcf_eval.eval_helper import chrom_evaluate_variants,ChromVariantStats,_genotype_concordance_dict
from smashbenchmarking.parsers.genome import Genome

class VariantsTestCase(unittest.TestCase):
    def testInit(self):
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr19   11   .       ACT     A       20      PASS    .       GT      1/1\n
chr19   15   .       G      A       20       PASS    .      GT      1/1\n
chr19   16   .       A      ATCG        20      PASS    .       GT      1/1\n
chr20   22   .       ATT       A         20      PASS    .       GT      0/1\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        newvcf = vcf.Reader(vcf_io)
        newvars = Variants(newvcf, MAX_INDEL_LEN)
        self.assertEqual(len(newvars.chroms),2)
        self.assertEqual(newvars.var_num(VARIANT_TYPE.SNP),1)
        self.assertEqual(newvars.var_num(VARIANT_TYPE.INDEL_DEL),2)
        self.assertEqual(newvars.var_num(VARIANT_TYPE.INDEL_INS),1)

    def testKnownFalsePositives(self):
        vcf_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr1    7       .       G       .       20      PASS    .       GT       0/0\n
"""
        vcf_io = StringIO.StringIO(vcf_str)
        newvcf = vcf.Reader(vcf_io)
        newvars = Variants(newvcf,MAX_INDEL_LEN,knownFP=True)
        chromvars = newvars.on_chrom('chr1')
        self.assertEqual(chromvars.all_locations,[7])
        self.assertEqual(chromvars.all_variants[7].alt[0],'NONE')


class VariantsHelpersTestCase(unittest.TestCase):
    def testAggregate(self):
        # build two ChromVariantStats objects
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
        gtdict = _genotype_concordance_dict() # leave empty for now
        cvs2 = chrom_evaluate_variants(true_vars,pred_vars,100,100,get_reference(),50)
        true_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr3   3   .       G     A       20      PASS    .       GT      1/1\n
chr3   5   .       C     T       20      PASS    .       GT      1/1\n
"""
        pred_str = """##fileformat=VCFv4.0\n
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001\n
chr3   3   .       G     A       20      PASS    .       GT      1/1\n
chr3   4   .       T     A       20      PASS    .       GT      1/1\n
chr3   7   .       G     C       20      PASS    .       GT      1/1\n
"""
        true_vars = vcf_to_ChromVariants(true_str,'chr3')
        pred_vars = vcf_to_ChromVariants(pred_str,'chr3')
        cvs3 = chrom_evaluate_variants(true_vars,pred_vars,100,100,get_reference(),50)
        #cvs5 = ChromVariantStats(true_vars,pred_vars,[31],[49,79],[52],_genotype_concordance_dict())
        aggregator,errors = _aggregate([cvs2,cvs3])
        # test some sums
        self.assertEqual(cvs2.num_true[VARIANT_TYPE.SNP],2)
        self.assertEqual(cvs3.num_true[VARIANT_TYPE.SNP],2)
        self.assertEqual(aggregator(VARIANT_TYPE.SNP)['num_true'],4)
        self.assertEqual(cvs2.num_tp[VARIANT_TYPE.SNP],1)
        self.assertEqual(cvs3.num_tp[VARIANT_TYPE.SNP],1)
        self.assertEqual(aggregator(VARIANT_TYPE.SNP)['good_predictions'],2)
        # TODO: test error function



if __name__ == "__main__":
    unittest.main()
# from vcfsmash import VariantComparator

# VCF_HEADER = """
# ##fileformat=VCFv4.0
# ##source=VCFWriter
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  NA12878
# """

# class SNPTestCases(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_single_perfect_match(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tT\tC\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = true_var 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': 1,
#             'small_ins_stats': 0,
#             'small_del_stats': 0,
#             'small_oth_stats': 0,
#             'large_ins_stats': 0,
#             'large_del_stats': 0,
#             'large_oth_stats': 0,
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value)
#             self.assertEqual(getattr(vc,field)['num_pred'], value)
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_negative(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tT\tC\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = VCF_HEADER 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (1,0),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_positive(self):
#         true_var = VCF_HEADER
#         predicted_var = (VCF_HEADER +
#             'chr1\t1061166\trs11807848\tT\tC\t20\tPASS\t.\tGT\t1/0')
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (0,1),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)
    
#     def test_mixed(self):
#         true_vars = VCF_HEADER
#         predicted_vars = VCF_HEADER
#         true_positives = [
#             'Chromosome\t474\t.\tC\tT\t96\tPASS\t.\tGT:DP\t1/1:23',
#             'Chromosome\t480\t.\tG\tA\t93\tPASS\t.\tGT:DP\t1/1:22',
#             'Chromosome\t393\t.\tT\tG\t87\tPASS\t.\tGT:DP\t1/1:20',
#             'Chromosome\t507\t.\tC\tT\t123\tPASS\t.\tGT:DP\t1/1:32',
#             'Chromosome\t542\t.\tG\tT\t171\tPASS\t.\tGT:DP\t1/1:48',
#             'Chromosome\t559\t.\tT\tC\t141\tPASS\t.\tGT:DP\t1/1:38',
#             'Chromosome\t588\t.\tG\tA\t57\tPASS\t.\tGT:DP\t1/1:10',
#             'Chromosome\t591\t.\tC\tT\t48\tPASS\t.\tGT:DP\t1/1:7',
#             'Chromosome\t600\t.\tG\tA\t78\tPASS\t.\tGT:DP\t1/1:17',
#             'Chromosome\t621\t.\tT\tC\t150\tPASS\t.\tGT:DP\t1/1:41',
#             'Chromosome\t660\t.\tT\tC\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t747\t.\tA\tG\t111\tPASS\t.\tGT:DP\t1/1:28'
#             ]
#         false_negatives = [
#             'Chromosome\t759\t.\tC\tT\t114\tPASS\t.\tGT:DP\t1/1:29',
#             'Chromosome\t774\t.\tT\tC\t111\tPASS\t.\tGT:DP\t1/1:28',
#             'Chromosome\t810\t.\tG\tT\t138\tPASS\t.\tGT:DP\t1/1:37',
#             'Chromosome\t831\t.\tC\tT\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t843\t.\tG\tA\t141\tPASS\t.\tGT:DP\t1/1:38',
#             'Chromosome\t882\t.\tT\tC\t126\tPASS\t.\tGT:DP\t1/1:33',
#             'Chromosome\t900\t.\tA\tT\t111\tPASS\t.\tGT:DP\t1/1:28',
#             ] 
#         false_positives = [
#             'Chromosome\t909\t.\tC\tT\t84\tPASS\t.\tGT:DP\t1/1:19',
#             'Chromosome\t966\t.\tT\tC\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t981\t.\tT\tG\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t1011\t.\tT\tC\t90\tPASS\t.\tGT:DP\t1/1:21',
#             'Chromosome\t1020\t.\tC\tT\t66\tPASS\t.\tGT:DP\t1/1:13'
#             ]
#         true_vars += '\n'.join(true_positives)
#         true_vars += '\n' + '\n'.join(false_negatives) 
        
#         predicted_vars += '\n'.join(true_positives)
#         predicted_vars += '\n' + '\n'.join(false_positives) 
        
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_vars)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_vars)
        
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         num_true = len(true_positives) + len(false_negatives)
#         num_predicted = len(true_positives) + len(false_positives)

#         expected = {
#             'snp_stats': (num_true,num_predicted),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])

#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

# class InsertTestCases(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_single_perfect_match(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tT\tTC\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = true_var 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': 0,
#             'small_ins_stats': 1,
#             'small_del_stats': 0,
#             'small_oth_stats': 0,
#             'large_ins_stats': 0,
#             'large_del_stats': 0,
#             'large_oth_stats': 0,
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value)
#             self.assertEqual(getattr(vc,field)['num_pred'], value)
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_negative(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tT\tTC\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = VCF_HEADER 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (1,0),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_positive(self):
#         true_var = VCF_HEADER
#         predicted_var = (VCF_HEADER +
#             'chr1\t1061166\trs11807848\tT\tTC\t20\tPASS\t.\tGT\t1/0')
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (0,1),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)
    
#     def test_mixed(self):
#         true_vars = VCF_HEADER
#         predicted_vars = VCF_HEADER
#         true_positives = [
#             'Chromosome\t474\t.\tC\tCT\t96\tPASS\t.\tGT:DP\t1/1:23',
#             'Chromosome\t480\t.\tG\tGC\t93\tPASS\t.\tGT:DP\t1/1:22',
#             'Chromosome\t393\t.\tT\tTG\t87\tPASS\t.\tGT:DP\t1/1:20',
#             'Chromosome\t507\t.\tC\tCTCTA\t123\tPASS\t.\tGT:DP\t1/1:32',
#             'Chromosome\t542\t.\tG\tGTGCT\t171\tPASS\t.\tGT:DP\t1/1:48',
#             'Chromosome\t559\t.\tT\tTC\t141\tPASS\t.\tGT:DP\t1/1:38',
#             # following line may break vcf spec 
#             'Chromosome\t588\t.\tG\tAG\t57\tPASS\t.\tGT:DP\t1/1:10',
#             'Chromosome\t591\t.\tC\tCT\t48\tPASS\t.\tGT:DP\t1/1:7',
#             'Chromosome\t600\t.\tG\tGCTATCCTTTAGAA\t78\tPASS\t.\tGT:DP\t1/1:17',
#             'Chromosome\t621\t.\tT\tTC\t150\tPASS\t.\tGT:DP\t1/1:41',
#             'Chromosome\t660\t.\tT\tTTTAGGAGCCAGAGACCCATGGAGATTATACACACACACAC\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t747\t.\tA\tAG\t111\tPASS\t.\tGT:DP\t1/1:28'
#             ]
#         false_negatives = [
#             'Chromosome\t759\t.\tC\tCT\t114\tPASS\t.\tGT:DP\t1/1:29',
#             'Chromosome\t774\t.\tT\tTC\t111\tPASS\t.\tGT:DP\t1/1:28',
#             'Chromosome\t810\t.\tG\tGTTCT\t138\tPASS\t.\tGT:DP\t1/1:37',
#             'Chromosome\t831\t.\tC\tCT\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t843\t.\tG\tGAGATATTAACCACCACACCACGATA\t141\tPASS\t.\tGT:DP\t1/1:38',
#             'Chromosome\t882\t.\tT\tTCCC\t126\tPASS\t.\tGT:DP\t1/1:33',
#             'Chromosome\t900\t.\tA\tATACAT\t111\tPASS\t.\tGT:DP\t1/1:28',
#             ] 
#         false_positives = [
#             'Chromosome\t909\t.\tC\tCTAT\t84\tPASS\t.\tGT:DP\t1/1:19',
#             'Chromosome\t966\t.\tT\tTACAGAC\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t981\t.\tT\tTATATATAG\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t1011\t.\tT\tTC\t90\tPASS\t.\tGT:DP\t1/1:21',
#             'Chromosome\t1020\t.\tC\tCT\t66\tPASS\t.\tGT:DP\t1/1:13'
#             ]
#         true_vars += '\n'.join(true_positives)
#         true_vars += '\n' + '\n'.join(false_negatives) 
        
#         predicted_vars += '\n'.join(true_positives)
#         predicted_vars += '\n' + '\n'.join(false_positives) 
        
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_vars)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_vars)
        
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         num_true = len(true_positives) + len(false_negatives)
#         num_predicted = len(true_positives) + len(false_positives)

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (num_true,num_predicted),
#             'small_del_stats': (0,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

# class InsertTestCases(unittest.TestCase):
#     def setUp(self):
#         pass

#     def tearDown(self):
#         pass

#     def test_single_perfect_match(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tTC\tT\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = true_var 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': 0,
#             'small_ins_stats': 0,
#             'small_del_stats': 1,
#             'small_oth_stats': 0,
#             'large_ins_stats': 0,
#             'large_del_stats': 0,
#             'large_oth_stats': 0,
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value)
#             self.assertEqual(getattr(vc,field)['num_pred'], value)
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_negative(self):
#         true_var = (VCF_HEADER + 
#             'chr1\t1061166\trs11807848\tTC\tT\t20\tPASS\t.\tGT\t1/0')
#         predicted_var = VCF_HEADER 
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (1,0),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

#     def test_single_false_positive(self):
#         true_var = VCF_HEADER
#         predicted_var = (VCF_HEADER +
#             'chr1\t1061166\trs11807848\tATTTATACATA\tA\t20\tPASS\t.\tGT\t1/0')
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_var)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_var)
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (0,1),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])
        
#         os.remove(true_vcf)
#         os.remove(predicted_vcf)
    
#     def test_mixed(self):
#         true_vars = VCF_HEADER
#         predicted_vars = VCF_HEADER
#         true_positives = [
#             'Chromosome\t474\t.\tCT\tC\t96\tPASS\t.\tGT:DP\t1/1:23',
#             'Chromosome\t480\t.\tGCTC\tG\t93\tPASS\t.\tGT:DP\t1/1:22',
#             'Chromosome\t393\t.\tTG\tT\t87\tPASS\t.\tGT:DP\t1/1:20',
#             'Chromosome\t507\t.\tCTACCATATA\tC\t123\tPASS\t.\tGT:DP\t1/1:32',
#             'Chromosome\t542\t.\tGA\tG\t171\tPASS\t.\tGT:DP\t1/1:48',
#             'Chromosome\t559\t.\tTAAAAAA\tT\t141\tPASS\t.\tGT:DP\t1/1:38',
#             'Chromosome\t588\t.\tAG\tA\t57\tPASS\t.\tGT:DP\t1/1:10',
#             'Chromosome\t591\t.\tCTATATC\tC\t48\tPASS\t.\tGT:DP\t1/1:7',
#             'Chromosome\t600\t.\tGTTTAATATAAAA\tG\t78\tPASS\t.\tGT:DP\t1/1:17',
#             'Chromosome\t621\t.\tTC\tT\t150\tPASS\t.\tGT:DP\t1/1:41',
#             'Chromosome\t660\t.\tTCC\tT\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t747\t.\tAG\tA\t111\tPASS\t.\tGT:DP\t1/1:28'
#             ]
#         false_negatives = [
#             'Chromosome\t759\t.\tCT\tC\t114\tPASS\t.\tGT:DP\t1/1:29',
#             'Chromosome\t774\t.\tTC\tT\t111\tPASS\t.\tGT:DP\t1/1:28',
#             'Chromosome\t810\t.\tGTCTAATA\tG\t138\tPASS\t.\tGT:DP\t1/1:37',
#             'Chromosome\t831\t.\tCT\tC\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t843\t.\tGAGA\tG\t141\tPASS\t.\tGT:DP\t1/1:38',
#             'Chromosome\t882\t.\tTATATACCA\tT\t126\tPASS\t.\tGT:DP\t1/1:33',
#             'Chromosome\t900\t.\tACATAC\tA\t111\tPASS\t.\tGT:DP\t1/1:28',
#             ] 
#         false_positives = [
#             'Chromosome\t909\t.\tCTAT\tC\t84\tPASS\t.\tGT:DP\t1/1:19',
#             'Chromosome\t966\t.\tTCA\tT\t135\tPASS\t.\tGT:DP\t1/1:36',
#             'Chromosome\t981\t.\tTCCCCCCC\tT\t159\tPASS\t.\tGT:DP\t1/1:44',
#             'Chromosome\t1011\t.\tTT\tT\t90\tPASS\t.\tGT:DP\t1/1:21',
#             'Chromosome\t1020\t.\tCCCC\tC\t66\tPASS\t.\tGT:DP\t1/1:13'
#             ]
#         true_vars += '\n'.join(true_positives)
#         true_vars += '\n' + '\n'.join(false_negatives) 
        
#         predicted_vars += '\n'.join(true_positives)
#         predicted_vars += '\n' + '\n'.join(false_positives) 
        
#         true_vcf = 'true_vars.vcf'
#         predicted_vcf = 'predicted_vars.vcf'
#         with open(true_vcf, 'w') as true_file:
#             true_file.write(true_vars)
#         with open(predicted_vcf, 'w') as predicted_file:
#             predicted_file.write(predicted_vars)
        
#         vc = VariantComparator(true_vcf, predicted_vcf)
#         vc.check_concordance()

#         num_true = len(true_positives) + len(false_negatives)
#         num_predicted = len(true_positives) + len(false_positives)

#         expected = {
#             'snp_stats': (0,0),
#             'small_ins_stats': (0,0),
#             'small_del_stats': (num_true,num_predicted),
#             'small_oth_stats': (0,0),
#             'large_ins_stats': (0,0),
#             'large_del_stats': (0,0),
#             'large_oth_stats': (0,0),
#             }
    
#         for field, value in expected.items():
#             self.assertEqual(getattr(vc,field)['num_true'], value[0])
#             self.assertEqual(getattr(vc,field)['num_pred'], value[1])

#         os.remove(true_vcf)
#         os.remove(predicted_vcf)

# if __name__ == '__main__':
#     unittest.main()

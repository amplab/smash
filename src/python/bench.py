#!/usr/bin/env python

#Copyright (c) 2013, Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Evaluate a predicted VCF against a "true" VCF.

Segregate by SNPs, indels, SVs.
Within indels and SVs, segregate by insertions, deletions, and "other".

Take error rate for validation data for SNPs, indels, and SVs.
WARNING: error rates are applied separately to insertions and deletions.
"""

from __future__ import division

from sys import argv
import vcf

from parsers.genome import Genome
from vcf_eval.variants import Variants,evaluate_variants,output_errors
from vcf_eval.chrom_variants import VARIANT_TYPE
from optparse import OptionParser

MAX_INDEL_LEN = 50
"""Largest size a variant can be while still being considered an indel.

Source: Alkan et al., "Genome structural variation discovery and genotyping".
"""


def nonzero_float(n):
  return float(n) if n != 0 else 1.0


def interval(lower, upper):
  """Format a confidence interval, given lower and upper bounds."""
  assert lower <= upper
  mid = (lower + upper) / 2
  radius = upper - mid
  return "%.1f +/- %.4f" % (mid * 100.0, radius * 100.0)


def bound_recall(tp, fn, e):
  """Bound recall, given numbers of TP, FN, and validation errors.

  This is a non-obvious theorem in the benchmarking paper.
  """
  num_true = tp + fn
  if not num_true:
    return 0, 0
  lower = (tp - e) / num_true
  upper = (tp + e) / num_true if e <= fn else tp / (num_true - e)
  return lower, upper

def bound_precision(tp, fp, e):
  """Bound precision, given numbers of TP, FP, and validation errors.

  This is a non-obvious theorem in the benchmarking paper.
  """
  p = tp + fp
  if not p or e > p:
    return 0, 0
  return (tp - e) / p, (tp + e) / p
  

def print_snp_results(num_true, num_pred, num_fp, num_fn, num_ib, num_ig, nrd, known_fp_prec, err):
  print "\n-----------"
  print "SNP Results"
  print "-----------"
  print "# True = %d; # Predicted = %d" % (num_true, num_pred)
  tp = num_ig
  assert tp + num_fn <= num_true
  assert tp + num_fp <= num_pred
  print "\t# precision =", interval(*bound_precision(tp, num_fp, err))
  if known_fp_var:
      print("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
  print "\t# recall =", interval(*bound_recall(tp, num_fn, err))
  print "\t# allele mismatch = %d" % num_ib
  print "\t# correct = %d" % num_ig
  print "\t# missed = %d" % num_fn
  print "\tpercent correct ignoring allele = %.1f" % (100 * (num_ig + num_ib) / nonzero_float(num_true))
  print "\tpercent correct = %.1f" % (100 * num_ig / nonzero_float(num_true))
  print "\tnon reference discrepancy = %1f" % (100*nrd)

def print_snp_stats(stats, err):
  print_snp_results(stats['num_true'], stats['num_pred'],
                    stats['false_positives'], stats['false_negatives'],
                    stats['intersect_bad'], stats['good_predictions'],ratio(stats['nrd_wrong'],stats['nrd_total']),
                    1-ratio(stats['known_fp_calls'],stats['known_fp']),err)

def ratio(a,b,sig=5):
    if b == 0:
        if a > 0:
            return 1.0
        return 0.0
    return float(int(10**sig*(float(a)/b)))/10**sig

def print_sv_results(var_type_str, num_true, num_pred, num_fp, num_fn, num_mm, num_gp, nrd, known_fp_prec, err):
  print "\n\n------------------------"
  print "%s Results" % var_type_str
  print "------------------------"
  print "# True = %d; # Predicted = %d" % (num_true, num_pred)
  print "\t# precision =", interval(*bound_precision(num_gp, num_fp, err))
  if known_fp_var:
        print("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
  print "\t# recall =", interval(*bound_recall(num_true - num_fn, num_fn, err))
  #print "\t# multiple matches = %d" % num_mm
  print "\t# correct = %d" % num_gp
  print "\t# missed = %d" % num_fn
  print "\t# false pos = %d" % num_fp
  if num_true > 0:
    print "\tpercent correct = %.1f" % (100 * num_gp / nonzero_float(num_true))
    print "\tnon reference discrepancy = %.1f" % (100*nrd)
  # The issue of multiple matches is empirically negligible.
  # TODO: assumed that this doesn't happen now; better way would be to assume that equidistant predicted var is wrong
  assert num_gp + num_fn <= num_true
  assert num_gp + num_fp <= num_pred

def print_sv_stats(description, stats, err):
  print_sv_results(description, stats['num_true'], stats['num_pred'],
                   stats['false_positives'], stats['false_negatives'],
                   #stats['mult_matches'],
                   0, stats['good_predictions'], ratio(stats['nrd_wrong'],stats['nrd_total']),
                   1-ratio(stats['known_fp_calls'],stats['known_fp']), err)
  
def print_sv_other_results(var_type_str, num_true, num_pred):
  print "\n\n------------------------"
  print "%s Statistics" % var_type_str
  print "------------------------"
  print "# True = %d; # Predicted = %d" % (num_true, num_pred) 


def _parse_args():
    usage = "usage: %prog [options] true.vcf pred.vcf (ref.fa)"
    parser = OptionParser(usage)
    parser.add_option("--fp",dest="knownFP",action="store",help="An optional VCF of known false-positives")
    parser.add_option("--sv_bp",dest="sv_eps",action="store",type=int,help="The maximum distance between SV breakpoints for them to be considered the same event",default=100)
    parser.add_option("--snp_err",dest="snp_err_rate",type=float,default=0.0,help="The error rate of SNPs in the ground truth data")
    parser.add_option("--indel_err",dest="indel_err_rate",type=float,default=0.0,help="The error rate of INDELs in the ground truth data")
    parser.add_option("--sv_err",dest="sv_err_rate",type=float,default=0.0,help="The error rate of SVs in the ground truth data")
    parser.add_option("-w","--rescue_window_size",dest="window",type=int,default=50,help="The size of the window for rescuing")
    parser.add_option("--err_vcf",dest="err_vcf",action="store",help="An optional output VCF to hold detected false-negatives and false-positives")
    options,args = parser.parse_args()
    if len(args) < 2:
      parser.error("Must specify input VCFs")
    return (options,args)

#------------
# MAIN SCRIPT
#------------


# Read command line args and parse VCF files
options,args = _parse_args()

# Read command line args and parse VCF files
true_vcf = vcf.Reader(open(args[0], 'r'))
pred_vcf = vcf.Reader(open(args[1], 'r'))
known_fp_vcf = vcf.Reader(open(options.knownFP,'r')) if options.knownFP else None

sv_eps, snp_err_rate, indel_err_rate, sv_err_rate = options.sv_eps, options.snp_err_rate, options.indel_err_rate, options.sv_err_rate

def stdAbbrev(ctig):
    return ctig.split()[0]

ref = Genome(args[2],stdAbbrev) if len(args) >= 3 else None
window = options.window if ref else None

true_var = Variants(true_vcf, MAX_INDEL_LEN)
pred_var = Variants(pred_vcf, MAX_INDEL_LEN)
known_fp_var = Variants(known_fp_vcf,MAX_INDEL_LEN,knownFP=True) if known_fp_vcf else None

# Estimated total number of errors in validation data for SNPs, indels and SVs.
snp_err = true_var.var_num(VARIANT_TYPE.SNP) * snp_err_rate
indel_err = (true_var.var_num(VARIANT_TYPE.INDEL_INS) + true_var.var_num(VARIANT_TYPE.INDEL_DEL) + true_var.var_num(VARIANT_TYPE.INDEL_OTH)) * indel_err_rate
sv_err = (true_var.var_num(VARIANT_TYPE.SV_INS) + true_var.var_num(VARIANT_TYPE.SV_DEL) + true_var.var_num(VARIANT_TYPE.SV_OTH)) * sv_err_rate

overall_statistics,errors = evaluate_variants(true_var,pred_var,sv_eps,sv_eps,ref,window,known_fp_var)
print_snp_stats(overall_statistics(VARIANT_TYPE.SNP), snp_err)

def print_sv(var_type, description):
  assert 'INDEL' in var_type or 'SV' in var_type
  err = indel_err if 'INDEL' in var_type else sv_err
  print_sv_stats(description, overall_statistics(var_type), err)

def print_oth(var_type, description):
  print_sv_other_results(description, true_var.var_num(var_type),
                         pred_var.var_num(var_type))

print_sv(VARIANT_TYPE.INDEL_DEL, 'INDEL DELETION')
print_sv(VARIANT_TYPE.INDEL_INS, 'INDEL INSERTION')
print_oth(VARIANT_TYPE.INDEL_OTH, 'INDEL OTHER')
print_sv(VARIANT_TYPE.SV_DEL, 'SV DELETION'),
print_sv(VARIANT_TYPE.SV_INS, 'SV INSERTION'),
print_oth(VARIANT_TYPE.SV_OTH, 'SV OTHER')

if options.err_vcf :
    output_errors(errors,ref.keys() if ref != None else None, open(options.err_vcf,'w'))

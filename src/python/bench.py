#!/usr/bin/env python

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

"""Evaluate a predicted VCF against a "true" VCF.

Segregate by SNPs, indels, SVs.
Within indels and SVs, segregate by insertions, deletions, and "other".

Take error rate for validation data for SNPs, indels, and SVs.
WARNING: error rates are applied separately to insertions and deletions.
"""

from __future__ import division, print_function

import os
import vcf
import sys
import argparse

from parsers.genome import Genome
from vcf_eval.variants import Variants,evaluate_variants,output_errors
from vcf_eval.chrom_variants import VARIANT_TYPE

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

def print_snp_results(num_true, num_pred, num_fp, num_fn, num_ib, num_ig, nrd,
        known_fp_prec, err, known_fp_vars=False, writer=print):
    writer("\n-----------")
    writer("SNP Results")
    writer("-----------")
    writer("# True = %d; # Predicted = %d" % (num_true, num_pred))
    tp = num_ig
    assert tp + num_fn <= num_true
    assert tp + num_fp <= num_pred
    writer("\t# precision =" + interval(*bound_precision(tp, num_fp, err)))
    if known_fp_vars:
        writer("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
    writer("\t# recall =" + interval(*bound_recall(tp, num_fn, err)))
    writer("\t# allele mismatch = %d" % num_ib)
    writer("\t# correct = %d" % num_ig)
    writer("\t# missed = %d" % num_fn)
    writer("\tpercent correct ignoring allele = %.1f" % (100 * (num_ig + num_ib) / nonzero_float(num_true)))
    writer("\tpercent correct = %.1f" % (100 * num_ig / nonzero_float(num_true)))
    writer("\tnon reference discrepancy = %1f" % (100*nrd))

def print_snp_stats(stats, err, known_fp_vars=False, writer=print):
    print_snp_results(stats['num_true'], stats['num_pred'],
                    stats['false_positives'], stats['false_negatives'],
                    stats['intersect_bad'], stats['good_predictions'],
                    ratio(stats['nrd_wrong'],stats['nrd_total']),
                    1-ratio(stats['known_fp_calls'],stats['known_fp']),
                    err, known_fp_vars, writer=writer)

def ratio(a,b,sig=5):
    if b == 0:
        if a > 0:
            return 1.0
        return 0.0
    return float(int(10**sig*(float(a)/b)))/10**sig

def print_sv_results(var_type_str, num_true, num_pred, num_fp, num_fn, num_mm,
        num_gp, nrd, known_fp_prec, err, known_fp_vars=False, writer=print):
    writer("\n\n------------------------")
    writer("%s Results" % var_type_str)
    writer("------------------------")
    writer("# True = %d; # Predicted = %d" % (num_true, num_pred))
    writer("\t# precision =", interval(*bound_precision(num_gp, num_fp, err)))
    if known_fp_vars:
          writer("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
    writer("\t# recall =", interval(*bound_recall(num_true - num_fn, num_fn, err)))
    #writer( "\t# multiple matches = %d" % num_mm)
    writer("\t# correct = %d" % num_gp)
    writer("\t# missed = %d" % num_fn)
    writer("\t# false pos = %d" % num_fp)
    if num_true > 0:
      writer("\tpercent correct = %.1f" % (100 * num_gp / nonzero_float(num_true)))
      writer("\tnon reference discrepancy = %.1f" % (100*nrd))
    # The issue of multiple matches is empirically negligible.
    # TODO: assumed that this doesn't happen now; better way would be to assume that equidistant predicted var is wrong
    assert num_gp + num_fn <= num_true
    assert num_gp + num_fp <= num_pred

def print_sv_stats(description, stats, err, writer=print):
    print_sv_results(description, stats['num_true'], stats['num_pred'],
                   stats['false_positives'], stats['false_negatives'],
                   #stats['mult_matches'],
                   0, stats['good_predictions'], ratio(stats['nrd_wrong'],
                   stats['nrd_total']), 1-ratio(stats['known_fp_calls'],
                   stats['known_fp']), err, writer=writer)

def print_other_results(var_type, num_true, num_pred, writer=print):
  writer("\n\n------------------------")
  writer("%s Statistics" % var_type)
  writer("------------------------")
  writer("# True = %d; # Predicted = %d" % (num_true, num_pred))

def parse_args(params):

    def is_valid_file(parser, arg):
        if not os.path.exists(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            return arg

    parser = argparse.ArgumentParser(description="""
        SMaSH benchmark toolkit for variant calling.
        See smash.cs.berkeley.edu for more information, including usage
        """)

    parser.add_argument('true_vcf', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('predicted_vcf', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('reference', type=lambda fn: is_valid_file(parser, fn),
            nargs='?')

    parser.add_argument("--fp",dest="knownFP",action="store",
            help="An optional VCF of known false-positives")
    parser.add_argument("--max_indel_len",dest="max_indel_len",action="store",
            type=int, default=MAX_INDEL_LEN,
            help="The maximum length of an indel or a \"small\" variation")
    parser.add_argument("--sv_bp",dest="sv_eps",action="store",type=int,
            default=100,
            help="""The maximum distance between SV breakpoints for them
                    to be considered the same event""")
    parser.add_argument("--snp_err",dest="snp_err_rate",type=float,default=0.0,
            help="The error rate of SNPs in the ground truth data")
    parser.add_argument("--indel_err",dest="indel_err_rate",type=float,
            default=0.0,
            help="The error rate of indels in the ground truth data")
    parser.add_argument("--sv_err",dest="sv_err_rate",type=float,default=0.0,
            help="The error rate of SVs in the ground truth data")
    parser.add_argument("-w","--rescue_window_size",dest="window",type=int,
            default=50,help="The size of the window for rescuing")
    parser.add_argument("--err_vcf",dest="err_vcf",action="store",
            help="""An optional output VCF to hold detected
            false-negatives and false-positives""")
    args = parser.parse_args(params)
    return args

class VariantComparator:
    def __init__(self, true_vars, pred_vars, max_indel_len=50, snp_err_rate=0,
            indel_err_rate=0, sv_err_rate = 0, sv_max_dist = 100):
        self.max_indel_len = max_indel_len
        self.snp_err_rate = snp_err_rate
        self.indel_err_rate = indel_err_rate
        self.sv_err_rate = sv_err_rate
        self.sv_max_dist = sv_max_dist
        self._true_vars = vcf.Reader(filename=true_vars)
        self._pred_vars = vcf.Reader(filename=pred_vars)
        self.true_vars = Variants(self._true_vars, self.max_indel_len)
        self.pred_vars = Variants(self._pred_vars, self.max_indel_len)
        self.known_fp_vars = None
        self.reference = None
        self.window = None

    def add_reference(self, reference, window):
        self.reference = Genome(reference, abbreviate= lambda c: c.split()[0])
        self.window = window

    def add_known_false_positives(self, known_fp_vars):
        self._known_fp_vars = vcf.Reader(filename=known_fp_vars)
        self.known_fp_vars = Variants(self._known_fp_vars, self.max_indel_len,
                knownFp=True)

    def check_concordance(self, err_output_file=None):
        # Estimated total number of errors in validation data for SNPs, indels and SVs.
        snp_err = self.true_vars.var_num(VARIANT_TYPE.SNP) * self.snp_err_rate
        indel_err = sum([
            self.true_vars.var_num(VARIANT_TYPE.INDEL_INS),
            self.true_vars.var_num(VARIANT_TYPE.INDEL_DEL),
            self.true_vars.var_num(VARIANT_TYPE.INDEL_OTH),
            ]) * self.indel_err_rate
        sv_err = sum([
            self.true_vars.var_num(VARIANT_TYPE.SV_INS),
            self.true_vars.var_num(VARIANT_TYPE.SV_DEL),
            self.true_vars.var_num(VARIANT_TYPE.SV_OTH)
            ]) * self.sv_err_rate
        self.stat_reporter, self.errors = evaluate_variants(self.true_vars,
                self.pred_vars, self.sv_max_dist, self.sv_max_dist,
                self.reference, self.window, self.known_fp_vars)
        self.snp_stats = self.stat_reporter(VARIANT_TYPE.SNP)
        self.small_ins_stats = self.stat_reporter(VARIANT_TYPE.INDEL_INS)
        self.small_del_stats = self.stat_reporter(VARIANT_TYPE.INDEL_DEL)
        self.small_oth_stats = (
                self.true_vars.var_num(VARIANT_TYPE.INDEL_OTH),
                self.pred_vars.var_num(VARIANT_TYPE.INDEL_OTH)
        )
        self.large_ins_stats = self.stat_reporter(VARIANT_TYPE.SV_INS)
        self.large_del_stats = self.stat_reporter(VARIANT_TYPE.SV_DEL)
        self.large_oth_stats = (
                self.true_vars.var_num(VARIANT_TYPE.SV_OTH),
                self.pred_vars.var_num(VARIANT_TYPE.SV_OTH)
        )
        if err_output_file:
            ref_keys = self.reference.keys() if self.reference else None
            with open(err_output_file, 'w') as err_file:
                output_errors(self.errors, ref_keys, err_file)

    def print_stats(self, writer=print):
        indel_err = self.indel_err_rate
        sv_err = self.sv_err_rate
        print_snp_stats(self.snp_stats, self.snp_err_rate, self.known_fp_vars,
                writer=writer)
        print_sv_stats('Small Insertions', self.small_ins_stats, indel_err,
                writer=writer)
        print_sv_stats('Small Deletions', self.small_del_stats, indel_err,
                writer=writer)
        print_other_results('Small Others', self.small_oth_stats[0],
                self.small_oth_stats[1], writer=writer)
        print_sv_stats('Large Insertions', self.large_ins_stats, sv_err,
                writer=writer)
        print_sv_stats('Large Deletions', self.large_del_stats, sv_err,
                writer=writer)
        print_other_results('Large Others', self.large_oth_stats[0],
                self.large_oth_stats[1], writer=writer)

def main(params):
    args = parse_args(params)
    var_comparator = VariantComparator(args.true_vcf, args.predicted_vcf,
            args.max_indel_len, args.snp_err_rate, args.indel_err_rate,
            args.sv_err_rate, args.sv_eps)
    if args.reference:
        var_comparator.add_reference(args.reference)
    if args.knownFP:
        var_comparator.add_known_false_positives(args.knownFP)
    var_comparator.check_concordance(args.err_vcf)
    var_comparator.print_stats()

if __name__ == '__main__':
    main(params=sys.argv[1:])

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
import csv
import sys
import argparse

from parsers.genome import Genome
from vcf_eval.variants import Variants,evaluate_variants,output_errors
from vcf_eval.chrom_variants import VARIANT_TYPE
from vcf_eval.callset_helper import MAX_INDEL_LEN
from normalize_vcf import NormalizeIterator

# this needs to move to another class
tsv_header = ['VariantType','#True','#Pred','Precision','Recall','TP','FP','FN']

def tsv_row(variant_name,stats,err):
    return [variant_name,
    stats['num_true'],
    stats['num_pred'],
    interval(*bound_precision(stats['good_predictions'],stats['false_positives'],err)),
    interval(*bound_recall(stats['good_predictions'],stats['false_negatives'],err)),
    stats['good_predictions'],
    stats['false_negatives'],
    stats['false_positives']
    ]

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

def print_snp_results(num_true, num_pred, num_fp, num_fn, num_ib, num_ig, nrd, known_fp_prec, err, known_fp_vars=False):
    print("\n-----------")
    print("SNP Results")
    print("-----------")
    print("# True = %d; # Predicted = %d" % (num_true, num_pred))
    tp = num_ig
    assert tp + num_fn <= num_true
    assert tp + num_fp <= num_pred
    print("\t# precision =", interval(*bound_precision(tp, num_fp, err)))
    if known_fp_vars:
        print("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
    print("\t# recall =", interval(*bound_recall(tp, num_fn, err)))
    print("\t# allele mismatch = %d" % num_ib)
    print("\t# correct = %d" % num_ig)
    print("\t# missed = %d" % num_fn)
    print("\tpercent correct ignoring allele = %.1f" % (100 * (num_ig + num_ib) / nonzero_float(num_true)))
    print("\tpercent correct = %.1f" % (100 * num_ig / nonzero_float(num_true)))
    print("\tnon reference discrepancy = %1f" % (100*nrd))

def print_snp_stats(stats, err, known_fp_vars=False):
    print_snp_results(stats['num_true'], stats['num_pred'],
                    stats['false_positives'], stats['false_negatives'],
                    stats['intersect_bad'], stats['good_predictions'],ratio(stats['nrd_wrong'],stats['nrd_total']),
                    1-ratio(stats['known_fp_calls'],stats['known_fp']),err, known_fp_vars)

def ratio(a,b,sig=5):
    if b == 0:
        if a > 0:
            return 1.0
        return 0.0
    return float(int(10**sig*(float(a)/b)))/10**sig

def print_sv_results(var_type_str, num_true, num_pred, num_fp, num_fn, num_mm, num_gp, nrd, known_fp_prec, err, known_fp_vars=False):
    print("\n\n------------------------")
    print("%s Results" % var_type_str)
    print("------------------------")
    print("# True = %d; # Predicted = %d" % (num_true, num_pred))
    print("\t# precision =", interval(*bound_precision(num_gp, num_fp, err)))
    if known_fp_vars:
        print("\t# precision (known FP) = %.1f" % (100*known_fp_prec) )
    print("\t# recall =", interval(*bound_recall(num_true - num_fn, num_fn, err)))
    #print "\t# multiple matches = %d" % num_mm
    print("\t# correct = %d" % num_gp)
    print("\t# missed = %d" % num_fn)
    print("\t# false pos = %d" % num_fp)
    if num_true > 0:
        print("\tpercent correct = %.1f" % (100 * num_gp / nonzero_float(num_true)))
        print("\tnon reference discrepancy = %.1f" % (100*nrd))
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
    print("\n\n------------------------")
    print("%s Statistics" % var_type_str)
    print("------------------------")
    print("# True = %d; # Predicted = %d" % (num_true, num_pred))



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
    parser.add_argument("--normalize",action="store_true",
            help="Optionally normalize variants before evaluating them; requires reference file")
    parser.add_argument("--output",action="store",
            help="Specify output type: plain text or TSV", default="text")
    args = parser.parse_args(params)
    return args

def get_snp_err(true_vars, snp_err_rate):
    return true_vars.var_num(VARIANT_TYPE.SNP) * snp_err_rate

def get_indel_err(true_vars, indel_err_rate):
    return sum([
        true_vars.var_num(VARIANT_TYPE.INDEL_INS),
        true_vars.var_num(VARIANT_TYPE.INDEL_DEL),
        true_vars.var_num(VARIANT_TYPE.INDEL_OTH),
        ]) * indel_err_rate

def get_sv_err(true_vars, sv_err_rate):
    return sum([
        true_vars.var_num(VARIANT_TYPE.SV_INS),
        true_vars.var_num(VARIANT_TYPE.SV_DEL),
        true_vars.var_num(VARIANT_TYPE.SV_OTH)
        ]) * sv_err_rate

def main(params):
    args = parse_args(params)
    if args.normalize and not args.reference:
        print("Normalization requires a reference file.",file=sys.stderr)

    if args.reference:
        ref = Genome(args.reference,abbreviate= lambda ctig: ctig.split()[0])
        window = args.window
    else:
        ref = None
        window = None

    with open(args.true_vcf) as f:
        true_vcf = vcf.Reader(f)
        if args.normalize:
            true_vcf = NormalizeIterator(ref,true_vcf)
        true_vars = Variants(true_vcf, MAX_INDEL_LEN)
    with open(args.predicted_vcf) as f:
        pred_vcf = vcf.Reader(f)
        if args.normalize:
            pred_vcf = NormalizeIterator(ref,pred_vcf)
        pred_vars = Variants(pred_vcf, MAX_INDEL_LEN)



    if args.knownFP:
        with open(args.knownFP) as f:
            known_fp_vcf = vcf.Reader(f)
            known_fp_vars = Variants(known_fp_vcf,
                    MAX_INDEL_LEN, knownFp=True)
    else:
        known_fp_vars = None

    # Estimated total number of errors in validation data for SNPs, indels and SVs.
    snp_err = get_snp_err(true_vars,args.snp_err_rate)

    indel_err = get_indel_err(true_vars,args.indel_err_rate)

    sv_err = get_sv_err(true_vars,args.sv_err_rate)

    sv_eps = args.sv_eps

    stat_reporter, errors = evaluate_variants(
        true_vars,
        pred_vars,
        sv_eps,
        sv_eps,
        ref,
        window,
        known_fp_vars
        )

    if args.output == "tsv":
        tsvwriter = csv.writer(sys.stdout, delimiter='\t')
        tsvwriter.writerow(tsv_header)
        tsvwriter.writerow(tsv_row("SNP",stat_reporter(VARIANT_TYPE.SNP),snp_err))
        tsvwriter.writerow(tsv_row("Indel Deletions",stat_reporter(VARIANT_TYPE.INDEL_DEL),indel_err))
        tsvwriter.writerow(tsv_row("Indel Insertions",stat_reporter(VARIANT_TYPE.INDEL_INS),indel_err))
        tsvwriter.writerow(tsv_row("Indel Other",stat_reporter(VARIANT_TYPE.INDEL_OTH),indel_err))
        tsvwriter.writerow(tsv_row("SV Deletions",stat_reporter(VARIANT_TYPE.SV_DEL),sv_err))
        tsvwriter.writerow(tsv_row("SV Insertions",stat_reporter(VARIANT_TYPE.SV_INS),sv_err))
        tsvwriter.writerow(tsv_row("SV Other",stat_reporter(VARIANT_TYPE.SV_OTH),sv_err))
    else:
        snp_stats = stat_reporter(VARIANT_TYPE.SNP)
        print_snp_stats(snp_stats, snp_err, known_fp_vars)
        def print_sv(var_type, description):
            assert 'INDEL' in var_type or 'SV' in var_type
            err = indel_err if 'INDEL' in var_type else sv_err
            print_sv_stats(description, stat_reporter(var_type), err)
    
        def print_oth(var_type, description):
            print_sv_other_results(description, stat_reporter(var_type)['num_true'], stat_reporter(var_type)['num_pred'])
    
        print_sv(VARIANT_TYPE.INDEL_DEL, 'INDEL DELETION')
        print_sv(VARIANT_TYPE.INDEL_INS, 'INDEL INSERTION')
        print_oth(VARIANT_TYPE.INDEL_OTH, 'INDEL OTHER')
        print_sv(VARIANT_TYPE.SV_DEL, 'SV DELETION'),
        print_sv(VARIANT_TYPE.SV_INS, 'SV INSERTION'),
        print_oth(VARIANT_TYPE.SV_OTH, 'SV OTHER')

    if args.err_vcf :
        output_errors(errors,ref.keys() if ref != None else None, open(args.err_vcf,'w'))

if __name__ == '__main__':
    main(params=sys.argv[1:])

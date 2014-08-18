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
import datetime
import json

from parsers.genome import Genome
from vcf_eval.variants import Variants,evaluate_variants,output_annotated_variants,evaluate_low_memory
from vcf_eval.chrom_variants import VARIANT_TYPE
from vcf_eval.callset_helper import MAX_INDEL_LEN
from normalize_vcf import normalize

# metadata
SMASHVERSION = "1.0"
date_run = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# this needs to move to another class
def get_tsv_header(knownFP=False):
    if knownFP:
        return ['VariantType','#True','#Pred','Precision', 'FP Precision','Recall','TP','FP','FN','NonReferenceDiscrepancy']
    else:
        return ['VariantType','#True','#Pred','Precision','Recall','TP','FP','FN','NonReferenceDiscrepancy']

def tsv_row(variant_name,stats,err_rate,knownFP=False,hideFP=False):
    err = err_rate * stats['num_true']
    if knownFP:
        return [variant_name,
        stats['num_true'],
        stats['num_pred'],
        interval(*bound_precision(stats['good_predictions'],stats['false_positives'],err)),
        interval(*bound_precision(stats['good_predictions'],stats['calls_at_known_fp'],err)),
        interval(*bound_recall(stats['good_predictions'],stats['false_negatives'],err)),
        stats['good_predictions'],
        stats['calls_at_known_fp'],
        stats['false_negatives'],
        get_nrd(stats)
        ]
    else:
        if hideFP:
            precision = "--"
            false_positives = "--"
        else:
            precision = interval(*bound_precision(stats['good_predictions'],stats['false_positives'],err))
            false_positives = stats['false_positives']
        return [variant_name,
        stats['num_true'],
        stats['num_pred'],
        precision,
        interval(*bound_recall(stats['good_predictions'],stats['false_negatives'],err)),
        stats['good_predictions'],
        false_positives,
        stats['false_negatives'],
        get_nrd(stats)
        ]

def json_dict(stats,err_rate,knownFP=False,hideFP=False):
    err = err_rate * stats['num_true']
    # should this be ordered dict?
    json_dict = {}
    if not knownFP:
        if hideFP:
            json_dict['precision'] = "--"
            json_dict['false_positives'] = "--"
        else:
            json_dict['precision'] = interval(*bound_precision(stats['good_predictions'],stats['false_positives'],err))
            json_dict['false_positives'] = stats['false_positives']
    json_dict['true'] = stats['num_true']
    json_dict['predicted'] = stats['num_pred']
    json_dict['true_positives'] = stats['good_predictions']
    json_dict['false_negatives'] = stats['false_negatives']
    json_dict['non_reference_discrepancy'] = get_nrd(stats)
    return json_dict

def get_nrd(stats):
    if stats['num_true'] > 0:
        return ratio(stats['nrd_wrong'],stats['nrd_total'])*100
    else:
        return 0.0

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

def print_snp_results(num_true, num_pred, num_fp, num_fn, num_ib, num_ig, nrd, known_fp_calls, err_rate, known_fp_vars=False,hideFP=False):
    err = num_true * err_rate
    print("\n-----------")
    print("SNP Results")
    print("-----------")
    print("# True = %d; # Predicted = %d" % (num_true, num_pred))
    tp = num_ig
    assert tp + num_fn <= num_true
    assert tp + num_fp <= num_pred
    if not hideFP:
        print("\t# precision =", interval(*bound_precision(tp, num_fp, err)))
    if known_fp_vars:
        print("\t# precision (known FP) =", interval(*bound_precision(tp,known_fp_calls,err)) )
    print("\t# recall =", interval(*bound_recall(tp, num_fn, err)))
    print("\t# allele mismatch = %d" % num_ib)
    print("\t# correct = %d" % num_ig)
    print("\t# missed = %d" % num_fn)
    print("\tpercent correct ignoring allele = %.1f" % (100 * (num_ig + num_ib) / nonzero_float(num_true)))
    print("\tpercent correct = %.1f" % (100 * num_ig / nonzero_float(num_true)))
    print("\tnon reference discrepancy = %1f" % (100*nrd))

def print_snp_stats(stats, err, known_fp_vars=False,hideFP=False):
    print_snp_results(stats['num_true'], stats['num_pred'],
                    stats['false_positives'], stats['false_negatives'],
                    stats['intersect_bad'], stats['good_predictions'],ratio(stats['nrd_wrong'],stats['nrd_total']),
                    stats['known_fp_calls'],err, known_fp_vars,hideFP)

def ratio(a,b,sig=5):
    if b == 0:
        if a > 0:
            return 1.0
        return 0.0
    return float(int(10**sig*(float(a)/b)))/10**sig

def print_sv_results(var_type_str, num_true, num_pred, num_fp, num_fn, num_mm, num_gp, nrd, known_fp_calls, err_rate, known_fp_vars=False,hideFP=False):
    err = err_rate * num_true
    print("\n\n------------------------")
    print("%s Results" % var_type_str)
    print("------------------------")
    print("# True = %d; # Predicted = %d" % (num_true, num_pred))
    if not hideFP:
        print("\t# precision =", interval(*bound_precision(num_gp, num_fp, err)))
    if known_fp_vars:
        print("\t# precision (known FP) = ", interval(*bound_precision(num_gp,known_fp_calls,err)) )
    print("\t# recall =", interval(*bound_recall(num_true - num_fn, num_fn, err)))
    #print "\t# multiple matches = %d" % num_mm
    print("\t# correct = %d" % num_gp)
    print("\t# missed = %d" % num_fn)
    if not hideFP:
        print("\t# false pos = %d" % num_fp)
    if num_true > 0:
        print("\tpercent correct = %.1f" % (100 * num_gp / nonzero_float(num_true)))
        print("\tnon reference discrepancy = %.1f" % (100*nrd))
    # The issue of multiple matches is empirically negligible.
    # TODO: assumed that this doesn't happen now; better way would be to assume that equidistant predicted var is wrong
    assert num_gp + num_fn <= num_true
    assert num_gp + num_fp <= num_pred

def print_sv_stats(description, stats, err,args):
    print_sv_results(description, stats['num_true'], stats['num_pred'],
                   stats['false_positives'], stats['false_negatives'],
                   #stats['mult_matches'],
                   0, stats['good_predictions'], ratio(stats['nrd_wrong'],stats['nrd_total']),
                   stats['known_fp_calls'], err,args.knownFP,args.hideFP)

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
    parser.add_argument('refindex', type=lambda fn: is_valid_file(parser,fn))

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
    parser.add_argument("--output_vcf",dest="output_vcf",action="store",
            help="""An optional output VCF to hold annotated variants from both source files""")
    parser.add_argument("--normalize",action="store_true",
            help="Optionally normalize variants before evaluating them; requires reference file")
    parser.add_argument("--output",action="store",
            help="Specify output type: plain text, TSV or JSON", default="text")
    parser.add_argument("--hide_fp",dest="hideFP",action="store_true",
            help="Don't show FP-related stats (for non-comprehensive ground truth files")
    args = parser.parse_args(params)
    return args

def get_text_header(params):
    return "# SMaSH version %s, run %s\n# cmdline args: %s" % (SMASHVERSION,date_run," ".join(params))

def add_json_header(params, output_dict):
    output_dict["SMaSH version"] = SMASHVERSION
    output_dict["Date run"] = date_run
    output_dict["cmdline args"] = " ".join(params)
    return output_dict

def get_vcf_header_lines(params):
    return ["##SMaSH version %s" % SMASHVERSION, "##Date run %s" % date_run, "##cmdline args: %s" % " ".join(params)]

def get_contig_lookup(fai):
    contig_list = []
    with open(fai,'r') as f:
        for line in f:
            contig_list.append(line.split('\t')[0])
    contig_lookup = {contig:index for index,contig in enumerate(contig_list)}
    max_index = max(contig_lookup.values())
    contig_lookup[None] = max_index + 100
    return contig_lookup

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

    contig_lookup = get_contig_lookup(args.refindex)

    true_vcf = vcf.Reader(open(args.true_vcf))
    if args.normalize:
        true_vcf = normalize(ref, true_vcf)
    pred_vcf = vcf.Reader(open(args.predicted_vcf))
    if args.normalize:
        pred_vcf = normalize(ref,pred_vcf)
    known_fp_vars = None # punting on known fp support for now

    sv_eps = args.sv_eps

    if args.output_vcf:
        outVCF = open(args.output_vcf,'w')
        outVCF.write("##fileformat=VCFv4.1\n")
        for s in get_vcf_header_lines(params):
            outVCF.write(s + "\n")
        outVCF.write("##INFO=<ID=smash_type,Type=String,Description=\"classify variant as TP,FP,FN,or rescued\">\n")
        outVCF.write("##INFO=<ID=source_file,Type=Integer,Description=\"variant originally in first or second vcf passed to SMaSH\">\n")
        outVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    else:
        outVCF = None

    stat_reporter = evaluate_low_memory(
        true_vcf,
        pred_vcf,
        sv_eps,
        sv_eps,
        ref,
        window,
        MAX_INDEL_LEN,
        contig_lookup,
        outVCF,
        known_fp_vars
        )

    if args.output == "tsv":
        print(get_text_header(params),file=sys.stdout)
        tsvwriter = csv.writer(sys.stdout, delimiter='\t')
        tsvwriter.writerow(get_tsv_header(args.knownFP))
        tsvwriter.writerow(tsv_row("SNP",stat_reporter(VARIANT_TYPE.SNP),args.snp_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("Indel Deletions",stat_reporter(VARIANT_TYPE.INDEL_DEL),args.indel_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("Indel Insertions",stat_reporter(VARIANT_TYPE.INDEL_INS),args.indel_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("Indel Inversions",stat_reporter(VARIANT_TYPE.INDEL_INV),args.indel_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("Indel Other",stat_reporter(VARIANT_TYPE.INDEL_OTH),args.indel_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("SV Deletions",stat_reporter(VARIANT_TYPE.SV_DEL),args.sv_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("SV Insertions",stat_reporter(VARIANT_TYPE.SV_INS),args.sv_err_rate,args.knownFP,args.hideFP))
        tsvwriter.writerow(tsv_row("SV Other",stat_reporter(VARIANT_TYPE.SV_OTH),args.sv_err_rate,args.knownFP,args.hideFP))
    elif args.output == "json":
        output_dict = add_json_header(params,{})
        output_dict["SNP"] = json_dict(stat_reporter(VARIANT_TYPE.SNP),args.snp_err_rate,args.knownFP,args.hideFP)
        output_dict["Indel Deletions"] = json_dict(stat_reporter(VARIANT_TYPE.INDEL_DEL),args.indel_err_rate,args.knownFP,args.hideFP)
        output_dict["Indel Insertions"] = json_dict(stat_reporter(VARIANT_TYPE.INDEL_INS),args.indel_err_rate,args.knownFP,args.hideFP)
        output_dict["Indel Inversions"] = json_dict(stat_reporter(VARIANT_TYPE.INDEL_INV),args.indel_err_rate,args.knownFP,args.hideFP)
        output_dict["Indel Other"] = json_dict(stat_reporter(VARIANT_TYPE.INDEL_OTH),args.indel_err_rate,args.knownFP,args.hideFP)
        output_dict["SV Deleitions"] = json_dict(stat_reporter(VARIANT_TYPE.SV_DEL),args.sv_err_rate,args.knownFP,args.hideFP)
        output_dict["SV Insertions"] = json_dict(stat_reporter(VARIANT_TYPE.SV_INS),args.sv_err_rate,args.knownFP,args.hideFP)
        json.dump(output_dict,sys.stdout)
    else:
        print(get_text_header(params),file=sys.stdout)
        snp_stats = stat_reporter(VARIANT_TYPE.SNP)
        print_snp_stats(snp_stats, args.snp_err_rate, known_fp_vars,args.hideFP)
        def print_sv(var_type, description,args):
            assert 'INDEL' in var_type or 'SV' in var_type
            err_rate = args.indel_err_rate if 'INDEL' in var_type else args.sv_err_rate
            print_sv_stats(description, stat_reporter(var_type), err_rate,args)
    
        def print_oth(var_type, description):
            print_sv_other_results(description, stat_reporter(var_type)['num_true'], stat_reporter(var_type)['num_pred'])
    
        print_sv(VARIANT_TYPE.INDEL_DEL, 'INDEL DELETION',args)
        print_sv(VARIANT_TYPE.INDEL_INS, 'INDEL INSERTION',args)
        print_sv(VARIANT_TYPE.INDEL_INV, 'INDEL INVERSION',args)
        print_oth(VARIANT_TYPE.INDEL_OTH, 'INDEL OTHER')
        print_sv(VARIANT_TYPE.SV_DEL, 'SV DELETION',args),
        print_sv(VARIANT_TYPE.SV_INS, 'SV INSERTION',args),
        print_oth(VARIANT_TYPE.SV_OTH, 'SV OTHER')

if __name__ == '__main__':
    main(params=sys.argv[1:])

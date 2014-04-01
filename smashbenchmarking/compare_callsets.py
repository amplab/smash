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

"""
Compare two callsets against each other where neither is known to be true
"""
# NB: since we're using methods for comparing a callset with ground truth
# some objects will have misleading names like "true positive", "false negative", etc.

# This is also an exercise in better organizing common code from bench.py, which is kind of a mess

from __future__ import division

import os
import vcf
import sys
import argparse

from parsers.genome import Genome
from vcf_eval.chrom_variants import VARIANT_TYPE
from vcf_eval.variants import Variants,evaluate_variants
from vcf_eval.callset_helper import MAX_INDEL_LEN

def perc(first,second):
    if second == 0:
        return 0.0
    else:
        return first/second*100

def print_results(stats,type_name):
    num_one = stats['num_true']
    num_two = stats['num_pred']
    print "\n-----------"
    print("{0} Results".format(type_name))
    print "-----------"
    print("# Callset 1 = {0}; # Callset 2 = {1}".format(num_one, num_two))
    intersection_totals = stats['good_predictions']
    print("# Intersection = {0}".format(intersection_totals))
    print("Callset 1: {0:.2f}% of calls in common with callset 2".format(perc(intersection_totals,num_one)))
    print("Callset 2: {0:.2f}% of calls in common with callset 1".format(perc(intersection_totals,num_two)))

# un-DRY; see similar function in bench.py
def parse_args(params):
    def is_valid_file(parser, arg):
        if not os.path.exists(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            return arg
    parser = argparse.ArgumentParser(description="""
        SMaSH benchmark toolkit for comparing two variant callsets.
        See smash.cs.berkeley.edu for more information, including usage
        """)
    parser.add_argument('callset1_vcf', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('callset2_vcf', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('reference', type=lambda fn: is_valid_file(parser, fn),
            nargs='?')
    parser.add_argument("--sv_bp",dest="sv_eps",action="store",type=int,
        default=100,
        help="""The maximum distance between SV breakpoints for them
                to be considered the same event""")
    parser.add_argument("-w","--rescue_window_size",dest="window",type=int,
        default=50,help="The size of the window for rescuing")
    args = parser.parse_args(params)
    return args

# also remarkably similar to the main in bench
def main(params):
    args = parse_args(params)
    with open(args.callset1_vcf) as f:
        callset1_vcf = vcf.Reader(f)
        one_vars = Variants(callset1_vcf,MAX_INDEL_LEN)
    with open(args.callset2_vcf) as f:
        callset2_vcf = vcf.Reader(f)
        two_vars = Variants(callset2_vcf,MAX_INDEL_LEN)
    if args.reference:
        ref = Genome(args.reference,abbreviate= lambda ctig: ctig.split()[0])
        window = args.window
    else:
        ref = None
        window = None

    stat_reporter, errors = evaluate_variants(
        one_vars,
        two_vars,
        args.sv_eps,
        args.sv_eps,
        ref,
        window
    )

    print_results(stat_reporter(VARIANT_TYPE.SNP),"SNP")
    print_results(stat_reporter(VARIANT_TYPE.INDEL_DEL),"INDEL DELETIONS")
    print_results(stat_reporter(VARIANT_TYPE.INDEL_INS),"INDEL INSERTIONS")
    print_results(stat_reporter(VARIANT_TYPE.INDEL_OTH),"INDEL OTHER")
    print_results(stat_reporter(VARIANT_TYPE.SV_DEL),"SV DELETIONS")
    print_results(stat_reporter(VARIANT_TYPE.SV_INS),"SV INSERTIONS")
    print_results(stat_reporter(VARIANT_TYPE.SV_OTH),"SV OTHER")


if __name__ == '__main__':
    main(params=sys.argv[1:])

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

"""Put a VCF file in a canonical form.

Three forms of canonicalization:
* Merge nearby variants.
* Chop off redundant bases on the right.
* Left-normalize.

Keep stats on how far left-normalizing.

TODO: preserve PASS versus not, annotation.

TODO: generalize for genomes other than mouse. (heterozygosity)

TODO: include merging in this script rather than a separate script.
"""


from __future__ import print_function

import sys,os
import argparse
import vcf
from itertools import groupby

from parsers.genome import Genome
from parsers.util import infoToStr,addInfoEntry,strToAlts
from vcf_eval.chrom_variants import is_sv

# TODO
# two iterators should become one
# better organization of these methods
# check asserts/exceptions to fail gracefully

info_norm_tag = "OP"
info_norm_header_line = """##INFO=<ID=""" + info_norm_tag + """,Number=1,Type=Integer,Description="Original position before normalization">"""

left_slides = []
same_position_discards = 0

def add_chrom_var(record, var_dict,ref_genome,verbose=False):
    # in case of position conflict, resolve
    while var_dict.has_key(record.POS):
        # if var in dict has slid and variant has not, switch them
        if not record.INFO.has_key(info_norm_tag) and var_dict[record.POS].INFO.has_key(info_norm_tag):
            var_temp = var_dict[record.POS]
            var_dict[record.POS] = record
            record = var_temp
        elif record.INFO.has_key(info_norm_tag):
            try:
                record = shift_until_not_overlapping(var_dict[record.POS],record,ref_genome)
            except AssertionError:
                print("failed denorm at %s %d" % (record.CHROM, record.POS), file=sys.stderr)
                return
        else:
            global same_position_discards 
            same_position_discards += 1
            if verbose:
                print("Variant already exists at %d; variant %d %s:%s will not be evaluated" % \
                    (record.POS,record.POS,record.REF,str(record.ALT)), file=sys.stderr)
            # same_position_discards.append(1)
            return # original VCF had variants at same position, discard all but first variant
    var_dict[record.POS] = record

def find_redundancy(strings):
    """Return the length of the longest common proper suffix.

    Assume there's at least two strings, and all are nonempty,
    and not all are equal.
    """
    min_len = min(map(len, strings))
    redundancy = 0
    while redundancy < min_len - 1:
        chars = map(lambda s: s[-(redundancy + 1)], strings)
        if len(set(chars)) == 1:
            redundancy += 1
        else:
            break
    return redundancy

def genotype(record):
    if not record.samples:
        return "."
    # return {0 : "0/0", 1 : "0/1", 2: "1/1", None : "."}[vcfrecord.samples[0].gt_type]
    # PyVCF's gt_type field only contains the values above. Pass the actual gt string through
    # to avoid converting other values, e.g. "2/1" to "0/1"
    if record.samples[0].gt_type == 1:
        return record.samples[0].gt_nums
    return {0 : "0/0", 1 : "0/1", 2: "1/1", None : "."}[record.samples[0].gt_type]

def mean(l):
    return float(sum(l))/len(l) if len(l) > 0 else float('nan')

def show_slides(slides):
    nonzero_slides = filter(lambda s: s != 0, slides)
    if len(slides) == 0:
        meanslides = 0
    else:
        meanslides = mean(slides)
    print("Average slide:", meanslides, file=sys.stderr)
    print("Fraction of variants that slide:",
          float(len(nonzero_slides)) / (1+len(slides)), file=sys.stderr)
    if nonzero_slides:
        print("Average nontrivial slides:", mean(nonzero_slides),
              file=sys.stderr)
    print("Variants discarded because of position collisions: %d" % (same_position_discards),file=sys.stderr)

def left_normalize(ref_genome, chrom, pos, ref_allele, alts):
    """Slide left until the last base of all alleles isn't the same.

    WARNING: Assume 0-based coordinates.
    """
    def same_last_base(alleles):
        last_bases = map(lambda allele: allele[-1], alleles)
        return len(set(last_bases)) == 1
    orig_pos = pos
    while same_last_base(alts + [ref_allele]):
        pos -= 1
        prev_base = ref_genome.ref(chrom, pos, pos + 1)
        def slide(allele):
            return prev_base + allele[:-1]
        ref_allele = slide(ref_allele)
        alts = map(slide, alts)
    left_slides.append(orig_pos - pos)
    return pos, ref_allele, alts

def shift_until_not_overlapping(var_one,var_two,ref_genome):
    # change to 0-based coord
    one_pos = var_one.POS - 1
    two_pos = var_two.POS - 1
    def same_last_base(alleles):
        last_bases = map(lambda allele: allele[-1], alleles)
        return len(set(last_bases)) == 1
    def same_first_base(alleles):
        first_bases = map(lambda allele: allele[0], alleles)
        return len(set(first_bases)) == 1
    def slide(pos,ref_allele,sliding_allele):
        right_endpoint = pos + len(ref_allele) - 1
        right_base = ref_genome.ref(var_two.CHROM,right_endpoint,right_endpoint+1)
        return sliding_allele[1:] + right_base
    two_ref_allele = var_two.REF
    two_alt_alleles = map(lambda a: str(a), var_two.ALT)

    while two_pos < (one_pos + len(var_one.REF)):
        two_pos += 1
        two_alt_alleles = map(lambda a: slide(two_pos,two_ref_allele,a),two_alt_alleles)
        two_ref_allele = slide(two_pos,two_ref_allele,two_ref_allele)
        assert(same_last_base([two_ref_allele] + two_alt_alleles))
        assert(same_first_base([two_ref_allele] + two_alt_alleles))
    var_two.POS = two_pos + 1 # back to 1-based coord
    var_two.REF = two_ref_allele
    var_two.ALT = map(lambda a: vcf.model._Substitution(a),two_alt_alleles)
    return var_two

# this should be the only method anything outside this file deals with
# returns an iterator over the normalized VCF records
def normalize(ref_genome, reader, maxIndelLen = 50, cleanOnly = False, verbose=False):
    norm_iter = NormalizeStepOneIterator(ref_genome,reader,maxIndelLen,cleanOnly)
    if cleanOnly:
        return norm_iter
    else:
        # resolve any chrom/pos collisions and output
        return normalize_generator(ref_genome,norm_iter,verbose)

def normalize_generator(genome,norm_iter,verbose=False):
    for key, chrom_group in groupby(norm_iter,lambda r: r.POS):
        var_dict = {}
        for record in chrom_group:
            add_chrom_var(record,var_dict,genome,verbose)
        for pos in sorted(var_dict):
            yield var_dict[pos]

def keep_variant(record,maxIndelLen=50):
    if ( record.FILTER != [] and record.FILTER != "." and record.FILTER != "PASS" and record.FILTER != None or genotype(record) == "0/0"):
        return False # record filtered
    if ( genotype(record) == "." and not is_sv(record, maxIndelLen)):
        return False # filter variants without genotype information unless they are structural variants
    if ( "N" in record.REF or any(map(lambda a: "N" in str(a), record.ALT))):
        return False # filter out ref/alts with N base
    return True

def normalize_variant(record,ref_genome,cleanOnly=False):
    orig_pos = record.POS
    pos = record.POS - 1 # left normalize assumes 0-based coord
    contig = record.CHROM
    record.REF = str(record.REF.upper())
    record.ALT = map(lambda a: str(a).upper(), record.ALT)
    if cleanOnly:
        record.ALT = strToAlts(record.ALT)
        return record
    ref = record.REF
    alts = record.ALT
    all_alleles = alts + [ref]
    redundancy = find_redundancy(all_alleles)
    if redundancy:
        def chop(allele):
            return allele[:-redundancy]
        ref = chop(ref)
        alts = map(chop,alts)
    pos,ref,alts = left_normalize(ref_genome,contig,pos,ref,alts)
    record.POS = pos + 1 # restore to 1-based coord
    record.REF = ref # string
    record.ALT = strToAlts(alts)
    if orig_pos != record.POS:
        record.INFO = addInfoEntry({},info_norm_tag,orig_pos)
    return record

class NormalizeStepOneIterator:
    def __init__(self,genome,vcfiter,max_indel_length=50,cleanOnly=False):
        self.genome = genome
        self.vcfiter = vcfiter
        self.max_indel_length = max_indel_length
        self.cleanOnly = cleanOnly

    def __iter__(self):
        return self

    def next(self):
        record = self.vcfiter.next()
        while (not keep_variant(record,self.max_indel_length)):
            record = self.vcfiter.next()
        return normalize_variant(record,self.genome,self.cleanOnly)

def parse_args(params):

    def is_valid_file(parser, arg):
        if not os.path.exists(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            return arg

    parser = argparse.ArgumentParser(description="""
        SMaSH benchmark tool for normalizing VCF files.
        See smash.cs.berkeley.edu for more information, including usage
        """)

    parser.add_argument('vcf', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('reference', type=lambda fn: is_valid_file(parser, fn))
    parser.add_argument('max_indel',default=50,nargs="?")
    parser.add_argument('--cleanonly',dest="cleanonly",default=False,action="store_true",help="""
        Filter and standardize variants without left-normalization.
        """)
    parser.add_argument('--verbose',dest="verbose",action="store_true",default=False,help="""
        Emit warnings when filtering by position collision, etc.
        """)
    args = parser.parse_args(params)
    return args

def main():
    args = parse_args(sys.argv[1:])
    vcf_reader = vcf.Reader(open(args.vcf))
    ref = args.reference
    max_indel_length = args.max_indel
    cleanOnly = args.cleanonly
    verbose = args.verbose
    # TODO: add INFO line for slide field in normalized header
    # if cleanOnly:
    #     vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
    # else:
    #     vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
    # using PyVCF's writer preserves all other info in the record
    vcf_writer = vcf.Writer(sys.stdout,vcf_reader)
    normiter = normalize(Genome(ref,lambda t: t.split()[0]), vcf_reader, max_indel_length, cleanOnly, verbose)
    map(lambda r: vcf_writer.write_record(r), normiter)
    if not cleanOnly:
        show_slides(left_slides)

if __name__ == '__main__':
    main()
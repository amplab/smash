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

import sys

from numpy import mean,histogram
from collections import OrderedDict

import vcf

import parsers.vcfwriter
from parsers.genome import Genome
from vcf_eval.chrom_variants import is_sv
#from vcf_eval.normalize import find_redundancy

info_norm_tag = "OP"

normalize_header = """##fileformat=VCFv4.0
##source=VCFWriter
##INFO=<ID=OP,Number=1,Type=Integer,Description="Original position before normalization">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  """

#TODO: probably move info methods to parsers.util?
def infoToStr(info):
    if len(info) == 0:
        return '.'
    else:
        return ";".join(map(lambda (k,v): str(k)+"="+str(v), info.iteritems()))

def addInfoEntry(info,key,value):
    if type(info) == dict:
        return OrderedDict({key:value})
    else:
        return info.update({key:value})

def add_var(variant, var_dict):
    if var_dict.has_key(variant.CHROM):
        if var_dict[variant.CHROM].has_key(variant.POS):
            var_dict[variant.CHROM][variant.POS].append(variant)
        else:
            var_dict[variant.CHROM][variant.POS] = [variant]
    else:
        var_dict[variant.CHROM] = {variant.POS:[variant]}

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

def genotype(vcfrecord):
    if not vcfrecord.samples:
        return "."
    # return {0 : "0/0", 1 : "0/1", 2: "1/1", None : "."}[vcfrecord.samples[0].gt_type]
    # PyVCF's gt_type field only contains the values above. Pass the actual gt string through
    # to avoid converting other values, e.g. "2/1" to "0/1"
    if vcfrecord.samples[0].gt_type == 1:
        return vcfrecord.samples[0].gt_nums
    return {0 : "0/0", 1 : "0/1", 2: "1/1", None : "."}[vcfrecord.samples[0].gt_type]

def write(record, writer):
    return writer.write_record(record.CHROM, record.POS, '.',
                               record.REF, ','.join(map(lambda a: str(a),record.ALT)), record.samples[0].gt_nums,infoToStr(record.INFO)) # TODO: more gtypes.


left_slides = []


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
def normalize(ref_genome, reader, maxIndelLen = 50, cleanOnly = False):
    norm_iter = NormalizeStepOneIterator(ref_genome,reader,maxIndelLen,cleanOnly)
    if cleanOnly:
        return norm_iter
    else:
        var_dict = {}
        # iterate over to build chrom:pos dict
        for record in norm_iter:
            add_var(record,var_dict)
        # resolve any chrom/pos collisions and output
        return normalize_generator(ref_genome,var_dict)

def keep_variant(record,maxIndelLen=50):
    if ( record.FILTER != [] and record.FILTER != "." and record.FILTER != "PASS" and record.FILTER != None or genotype(record) == "0/0"):
        return False # record filtered
    if ( genotype(record) == "." and not is_sv(record, maxIndelLen)):
        return False # filter variants without genotype information unless they are structural variants
    return True

def normalize_variant(record,ref_genome,cleanOnly=False):
    orig_pos = record.POS
    pos = record.POS - 1 # left normalize assumes 0-based coord
    contig = record.CHROM
    record.REF = str(record.REF.upper())
    record.ALT = map(lambda a: str(a).upper(), record.ALT)
    if cleanOnly:
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
    record.ALT = map(lambda a: vcf.model._Substitution(a),alts) # for some reason alts are not strings in PyVCF?
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


def normalize_generator(genome,var_dict):
    for chrom in sorted(var_dict):
        for pos in sorted(var_dict[chrom]):
            if len(var_dict[chrom][pos]) == 1:
                yield var_dict[chrom][pos][0]
            else:
                resolved_vars = handle_position_collisions(var_dict[chrom][pos],genome)
                for v in resolved_vars:
                    yield v

def handle_position_collisions(variants,ref_genome):
    normed = []
    notnormed = []
    final_vars = []
    for r in variants:
        if r.INFO.has_key(info_norm_tag):
            normed.append(r)
        else:
            notnormed.append(r)
    # choose one of the not normed and throw the rest away
    if len(notnormed) != 0:
        basevar = notnormed[0]
        for v in notnormed[1:]:
            print("Variant already exists on %s at %d; discarding variant %s %d %s/%s" % (v.CHROM, v.POS, v.CHROM, v.POS, v.REF, str(v.ALT)),file=sys.stderr)
    # if no not normed, arbitrarily choose one to keep in initial position
    else:
        basevar = normed[0]
        normed = normed[1:]
    # then, place normed variants back after initial variant
    for r in normed:
        try:
            shiftedvar = shift_until_not_overlapping(basevar,r,ref_genome)
            final_vars.append(basevar)
            basevar = shiftedvar
        except AssertionError:
            print("failed denorm at %s %d" % (normed[0].CHROM, normed[0].POS), file=sys.stderr)
    final_vars.append(basevar)
    return final_vars

def main():
    vcf_reader = vcf.Reader(sys.stdin)
    ref = sys.argv[1]
    person = sys.argv[2]
    max_indel_length = 50
    if len(sys.argv) > 3:
        max_indel_length = sys.argv[3]
    cleanOnly = len(sys.argv) > 4 and sys.argv[4] == "cleanonly"
    if cleanOnly:
        vcf_writer = parsers.vcfwriter.VCFWriter(ref, person, sys.stdout)
    else:
        vcf_writer = parsers.vcfwriter.VCFWriter(ref,person,sys.stdout,normalize_header)
    normiter = normalize(Genome(ref,lambda t: t.split()[0]), vcf_reader, max_indel_length, cleanOnly)
    map(lambda r: write(r,vcf_writer), normiter)
    if not cleanOnly:
        show_slides(left_slides)


if __name__ == '__main__':
    main()

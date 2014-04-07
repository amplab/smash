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

from numpy import mean

import vcf

import parsers.vcfwriter
from parsers.genome import Genome
from vcf_eval.chrom_variants import is_sv

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

def write(record, pos, ref, alts, writer):
    return writer.write_record(record.CHROM, pos, '.',
                               ref, ','.join(alts), genotype(record)) # TODO: more gtypes.


left_slides = []


def show_slides(slides):
    nonzero_slides = filter(lambda s: s != 0, slides)
    print("Average slide:", mean(slides), file=sys.stderr)
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


def normalize(ref_genome, reader, writer, maxIndelLen = 50, cleanOnly = False):
    for record in reader:
        if ( record.FILTER != [] and record.FILTER != "." and record.FILTER != "PASS" and record.FILTER != None or genotype(record) == "0/0"):
            continue # record filtered
        if ( genotype(record) == "." and not is_sv(record, maxIndelLen)):
            continue # filter variants without genotype information unless they are structural variants
        pos = record.POS - 1
        contig = record.CHROM
        ref = str(record.REF.upper())
        alts = map(lambda a: str(a).upper(), record.ALT)
        all_alleles = alts + [ref]
        redundancy = find_redundancy(all_alleles)
        if cleanOnly:
            write(record,pos,ref,alts,writer)
            continue

        if redundancy:
            def chop(allele):
                return allele[:-redundancy]
            ref = chop(ref)
            alts = map(chop, alts)
        pos, ref, alts = left_normalize(ref_genome, record.CHROM,
                                        pos, ref, alts)
        write(record, pos, ref, alts, writer)
        prev_end = pos + len(ref)
        prev_contig = contig


def main():
    # ref_genome = genome.Genome(sys.argv[1])
    vcf_reader = vcf.Reader(sys.stdin)
    ref = sys.argv[1]
    person = sys.argv[2]
    max_indel_length = 50
    if len(sys.argv[3]) > 3:
        max_indel_length = sys.argv[3]
    cleanOnly = len(sys.argv) > 4 and sys.argv[4] == "cleanonly"
    vcf_writer = parsers.vcfwriter.VCFWriter(ref, person, sys.stdout)
    normalize(Genome(ref,lambda t: t.split()[0]), vcf_reader, vcf_writer, max_indel_length, cleanOnly)
    show_slides(left_slides)


if __name__ == '__main__':
    main()

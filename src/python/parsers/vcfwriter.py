"""Write a VCF file.

Intended application: benchmarking.
"""


from __future__ import print_function

import genome
import util


_anon_header = """##fileformat=VCFv4.0
##source=VCFWriter
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"""


def _must_prepend(ref, alts):
    """Return whether the preceding reference base must be prepended.

    Helper for satisfying the VCF spec.
    """
    alleles = alts + [ref]
    if any(not allele for allele in alleles):
        return True
    snp = all(len(allele) == 1 for allele in alleles)
    return not snp and any(allele[0] != ref[0] for allele in alts)


class VCFWriter:

    """VCF writer for a particular person and reference genome.

    Coordinates are "space-counted, zero-start" unless otherwise specified.
    <http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html>
    (Note that the VCF format itself is "base-counted, one-start".)
    """
    
    def __init__(self, reference_fasta, person, output):
        """Given a reference, a person, and a file-like object for output.

        Ideally the reference would be encoded in the FASTA.
        """
        self._ref_genome = genome.Genome(reference_fasta)
        self._output = output
        print(_anon_header + person, file=self._output)
        self._person = person

    def write_record(self, CHROM, POS, ID, REF, ALT, gtype):
        """Write a fully specified VCF record.

        WARNING: 'REF' isn't checked against the reference genome.

        Do nothing if 'REF' contains characters outside ACGT (notably N).
        Return whether anything was written.
        """
        write = util.is_proper_strand(REF)
        if write:
            QUAL = 20               # Default 1/100 error probability.
            FILTER = 'PASS'
            FORMAT = 'GT'
            print(CHROM, POS + 1, ID, REF, ALT, QUAL, FILTER, '.', FORMAT, gtype,
                  sep='\t', file=self._output)
        return write

    def write_deletion(self, CHROM, start, end, ID, gtype):
        """Write a deletion by looking up the deleted reference bases."""
        return self.write_deletion_with_insertion(CHROM, start, end,
                                                  ID, [''], gtype)
        
    def write_insertion(self, CHROM, start, inserted_sequence, ID, gtype):
        """Write an insertion by taking the preceding base as ref allele."""
        return self.write_deletion_with_insertion(CHROM, start, start,
                                                  ID, inserted_sequence, gtype)

    def write_deletion_with_insertion(self, CHROM, start, end, ID, alts,
                                      gtype):
        """Replace arbitrary sequence with arbitrary sequence.

        'alts' is a list of alleles.

        To conform to the VCF spec, prepend the preceding reference base to the
        REF and ALT alleles if any of them are empty, or if this variant isn't
        a SNP and not all alleles share a common first base.
        
        WARNING: If prepending, assume the base preceding the deletion wasn't
        in the reference allele of the previous variant, and that the deletion
        isn't at the start of a chromosome.
        """
        assert(0 <= start <= end)
        REF = self._ref_genome.ref(CHROM, start, end)
        def write(pos, ref, alts):
            ALT = ','.join(alts) if alts else '.'
            return self.write_record(CHROM, pos, ID, ref, ALT, gtype)
        if _must_prepend(REF, alts):
            assert(start)
            anchor = self._ref_genome.ref(CHROM, start - 1, start)
            return write(start - 1, anchor + REF, map(anchor.__add__, alts))
        else:
            return write(start, REF, alts)

    def write_alleles(self, CHROM, start, end, ID, alleles, phased=True):
        """Like 'write_deletion_with_insertion' but with a pair of alleles."""
        assert(1 <= len(alleles) <= 2)
        REF = self._ref_genome.ref(CHROM, start, end)
        distinct_alleles = [REF]
        for allele in alleles:
            if allele not in distinct_alleles:
                distinct_alleles.append(allele)
        allele_indices = map(distinct_alleles.index, alleles)
        sep = '|' if phased else '/'
        gtype = sep.join(map(str, allele_indices))
        alts = distinct_alleles[1:]
        return self.write_deletion_with_insertion(CHROM, start, end, ID, alts,
                                                  gtype)

    def write_inversion(self, CHROM, start, end, ID, gtype):
        assert(0 <= start < end)
        REF = self._ref_genome.ref(CHROM, start, end)
        ALT = REF[::-1]
        return self.write_record(CHROM, start, ID, REF, ALT, gtype)

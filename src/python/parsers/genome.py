import random

from pyfasta import fasta

import util


def _make_chrom_unabbrev(f, abbreviate):
    """Return a dictionary for unabbreviating chromosome names.

    Assume 'f' is a 'fasta.Fasta'.
    """
    chrom_unabbrev = {}
    for chrom in f.iterkeys():
        key = abbreviate(chrom)
        if key not in chrom_unabbrev:
            chrom_unabbrev[key] = chrom
        else:
            raise Exception('FASTA chromosome abbreviation collision.')
    return chrom_unabbrev


class Genome:

    """FASTA wrapper.

    <https://github.com/brentp/pyfasta>
    """

    def __init__(self, name, abbreviate=None):
        """Wrap pyfasta's 'Fasta' class.

        Optionally introduce a chromosome abbreviation function, 'abbreviate'.
        This feature is present in pyfasta's 'Fasta' class,
        but we reimplement it because it didn't seem to work properly,
        and we want to provide unabbreviation.
        """
        self._fasta = fasta.Fasta(name)
        self._abbrev = abbreviate
        if abbreviate:
            self._unabbrev = _make_chrom_unabbrev(self._fasta, abbreviate)

    def keys(self):
        return map(self._abbrev,self._fasta.keys())
    
    def unabbreviate(self, chrom):
        return self._unabbrev[chrom] if self._abbrev else chrom

    def _seq(self, chrom):
        return self._fasta[self.unabbreviate(chrom)]

    def chrom_length(self, chrom):
        return len(self._seq(chrom))

    def ref(self, chrom, start, end, orientation=True):
        """Return a 0-based slice of a reference chromosome.

        If 'orientation' is false, effectively reverse-complement the
        underlying sequence before slicing.
        """
        seq = self._seq(chrom)
        length = len(seq)
        if not orientation:
            start, end = length - end, length - start
        subseq = seq[start:end].upper()
        rc = util.reverse_complement_strand
        return subseq if orientation else rc(subseq)

    def random_sequence(self, chrom, length):
        """Return a random segment of given 'length' from given 'chrom'.

        Not quite random: the sequence won't have any "N".
        """
        while True:
            start = random.randint(0, len(self._fasta[chrom]) - length)
            seq = self.ref(chrom, start, start + length)
            if util.is_proper_strand(seq):
                return seq

    def close(self):
        pass


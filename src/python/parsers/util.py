"""Helpers.

TODO: Reduce code duplication with 'dbsnp/util.py'.
"""


import csv
import random


_BASES = frozenset(['A', 'T', 'C', 'G'])


def forall(pred, collection):
   for item in collection:
      if not pred(item):
         return False
   return True


def exists(pred, collection):
   for item in collection:
      if pred(item):
         return True
   return False


def is_proper_strand(strand):
   """Return whether a string consists of one or more DNA bases."""
   return strand and forall(_BASES.__contains__, strand)


def _is_strand(strand):
   return is_proper_strand(strand) or strand == '-'


def _is_alleles(alleles):
   strands = alleles.split('/')
   return strands and forall(_is_strand, strands)


_PAIRING = {'-':'-', 'A':'T', 'T':'A', 'C':'G', 'G':'C'}


def reverse_complement_strand(strand):
    rev_comp = []
    for base in strand:
       rev_comp.append(_PAIRING[base])
    rev_comp.reverse()
    return ''.join(rev_comp)


def reverse_complement(alleles):
   return '/'.join(map(_reverse_complement_strand, alleles.split('/')))


class MissingColumnError(Exception):
   pass


def parse_csv(filename, parse_row):
    """Process a CSV file by looking up entries by column name.

    For each row, 'parse_row' is called with single argument a lookup function
    for that row.  Multiple arguments can be given to the lookup function; the
    first argument corresponding to an existing column is used.
    """
    reader = csv.reader(open(filename, 'rb'), delimiter='\t')
    column_names = reader.next()
    for row in reader:
        def get_entry(*candidate_column_names):
           for column_name in candidate_column_names:
              if column_name in column_names:
                 return row[column_names.index(column_name)]
           raise MissingColumnError()
        parse_row(get_entry)

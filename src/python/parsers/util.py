#Copyright (c) 2013, Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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

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


"""Helper class for VCF evaluation."""

from eval_helper import chrom_evaluate_variants
from chrom_variants import ChromVariants
from itertools import chain

from collections import defaultdict

class Variants:

  """Variation across a genome, organized by chromosome."""
  
  def __init__(self, new_vcf, max_indel_len, **kwargs):
    self._args = kwargs
    self._max_indel_len = max_indel_len
    self._variation = {}
    for record in new_vcf:
      self._add_record(record)
    self._chroms = frozenset(self._variation.iterkeys())
    for c in self.chroms:
      self._variation[c].validate()

  @property
  def chroms(self):
    """Return the set of all chromosomes that have variants."""
    return self._chroms

  def _virgin_chrom(self, chrom):
    """Return an empty set of variants on a given chromosome."""
    return ChromVariants(chrom, self._max_indel_len, **self._args)

  def _add_record(self, record):
    chrom = record.CHROM
    if chrom not in self._variation:
      self._variation[chrom] = self._virgin_chrom(chrom)
    self._variation[chrom].add_record(record)
    
  def on_chrom(self, chrom):
    """Return all variants on 'chrom'.  'defaultdict'-like behavior."""
    return self._variation.get(chrom, self._virgin_chrom(chrom))
  
  def var_num(self, var_type):
    """Return the total number of variants of a given type."""
    def chrom_var_num(variants):
       return variants.var_num(var_type)
    return sum(map(chrom_var_num, self._variation.itervalues()))


def _aggregate(stats_by_chrom):
  """Combine dictionaries by summing their values."""
  # a bit tricky
  # return a function that takes a type and aggregates
  # the statistics by chromosome
  def aggregator(vartype):

   aggregate = defaultdict(int)
   for cv_stats in stats_by_chrom:
     stats = cv_stats.to_dict(vartype)
     for attribute in stats:
      aggregate[attribute] += stats[attribute]
   return dict(aggregate)

  def errors(chrom_order):
   if chrom_order == None:
       chrom_order = sorted(map(lambda t: t.chrom, stats_by_chrom))
   stats_by_chrom_dict = dict(map(lambda t: [t.chrom,t],stats_by_chrom) )
   var_iterator = []
   for chrom in filter(lambda t: t in stats_by_chrom_dict.keys(),chrom_order):
     cv_stats = stats_by_chrom_dict[chrom]
     fp_iter = cv_stats.false_positives._vcf_iterator("err_type=FP",cv_stats.false_positives.all_locations)
     fn_iter = cv_stats.false_negatives._vcf_iterator("err_type=FN",cv_stats.false_negatives.all_locations)
     if ( cv_stats.known_fp_variants != None ):
      known_fp_iter = cv_stats.known_fp_variants._vcf_iterator("err_type=FP_known",cv_stats.known_fp_variants.all_locations)
     else:
      known_fp_iter = [].__iter__()
     var_iterator = chain(var_iterator,vcf_by_position([fp_iter,fn_iter,known_fp_iter]))
   return var_iterator

  return aggregator,errors


def _eval_aggregate(true_variants, pred_variants, known_fp, evaluate):
  """Evaluate by chromosome and aggregate the results."""
  dicts = []
  for chrom in true_variants.chroms.union(pred_variants.chroms):
    chrom_true = true_variants.on_chrom(chrom)
    chrom_pred = pred_variants.on_chrom(chrom)
    chrom_known_fp = known_fp.on_chrom(chrom) if known_fp else None
    dicts.append(evaluate(chrom_true, chrom_pred,chrom_known_fp))
  return _aggregate(dicts)  

def evaluate_variants(true_variants,pred_variants,eps,eps_bp,ref,window,known_fp):
    def evaluate(true,pred,known_fp):
      return chrom_evaluate_variants(true,pred,known_fp,eps,eps_bp,ref,window)
    return _eval_aggregate(true_variants,pred_variants,known_fp,evaluate)

def output_errors(err_aggregate,contig_ordering,outVCF):
    # write the header
    outVCF.write("##fileformat=VCFv4.1\n")
    outVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for rec in err_aggregate(contig_ordering):
      outVCF.write("%s\n" % rec )

import heapq
def vcf_by_position(rec_iters):
     '''
     Given several iterators over VCFs (lines, not PyVCF records), return an iterator that presents them in sorted order
     '''
     def nextLoc(rec,default=None):
      p = int(rec.split("\t",2)[1])
      return (p,rec)
     pos_iters = map(lambda rec_iter: map(nextLoc,rec_iter).__iter__(),rec_iters)
     merge_iter = heapq.merge(*pos_iters)
     return map(lambda t: t[1], merge_iter).__iter__()

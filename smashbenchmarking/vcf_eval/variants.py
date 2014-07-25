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

"""Helper class for VCF evaluation."""

from eval_helper import chrom_evaluate_variants,chrom_var_stats_dict
from chrom_variants import ChromVariants,VARIANT_TYPE
from itertools import chain,groupby

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
      self._variation[c]._ensure_sorted()
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

# for stats: instead of returning a function and a whole-genome dict
# eval the function now to hold all stats in ONE dict, not broken out by chrom
# for vcfs likewise, get iter for 1 chrom and output, then release
# eventually these two methods completely replace the _aggregate method below

def aggregate_stats(genome_stats, chrom_stats):
  for vartype in VARIANT_TYPE:
    for attribute in chrom_stats.to_dict(vartype).keys():
      genome_stats[vartype][attribute] += chrom_stats.to_dict(vartype)[attribute]
  return genome_stats

def write_annotated_var(writer,cv_stats):
  var_iterator = []
  tp_iter = cv_stats.true_positives._vcf_iterator("source_file=1;smash_type=TP",cv_stats.true_positives.all_locations)
  fp_iter = cv_stats.false_positives._vcf_iterator("source_file=2;smash_type=FP",cv_stats.false_positives.all_locations)
  fn_iter = cv_stats.false_negatives._vcf_iterator("source_file=1;smash_type=FN",cv_stats.false_negatives.all_locations)
  rescued_iter = cv_stats.rescued_vars._vcf_iterator("source_file=2;smash_type=rescued",cv_stats.rescued_vars.all_locations)
  if (cv_stats.known_fp_variants is not None):
    # does this annotation make sense?
    known_fp_iter = cv_stats.known_fp_variants._vcf_iterator("err_type=FP_known",cv_stats.known_fp_variants.all_locations)
  else:
    known_fp_iter = [].__iter__()
  var_iterator = sorted(chain(var_iterator,vcf_by_position([tp_iter,fp_iter,fn_iter,known_fp_iter,rescued_iter])),key=lambda r: int(r.split("\t",2)[1]))
  map(lambda r: writer.write("%s\n" % r),var_iterator)

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

  def annotated_var_iter(chrom_order):
   if chrom_order == None:
       chrom_order = sorted(map(lambda t: t.chrom, stats_by_chrom))
   stats_by_chrom_dict = dict(map(lambda t: [t.chrom,t],stats_by_chrom) )
   var_iterator = []
   for chrom in filter(lambda t: t in stats_by_chrom_dict.keys(),chrom_order):
     cv_stats = stats_by_chrom_dict[chrom]
     tp_iter = cv_stats.true_positives._vcf_iterator("source_file=1;smash_type=TP",cv_stats.true_positives.all_locations)
     fp_iter = cv_stats.false_positives._vcf_iterator("source_file=2;smash_type=FP",cv_stats.false_positives.all_locations)
     fn_iter = cv_stats.false_negatives._vcf_iterator("source_file=1;smash_type=FN",cv_stats.false_negatives.all_locations)
     rescued_iter = cv_stats.rescued_vars._vcf_iterator("source_file=2;smash_type=rescued",cv_stats.rescued_vars.all_locations)
     if ( cv_stats.known_fp_variants != None ):
      known_fp_iter = cv_stats.known_fp_variants._vcf_iterator("err_type=FP_known",cv_stats.known_fp_variants.all_locations)
     else:
      known_fp_iter = [].__iter__()
     var_iterator = sorted(chain(var_iterator,vcf_by_position([tp_iter,fp_iter,fn_iter,known_fp_iter,rescued_iter])),key=lambda rec: int(rec.split("\t",2)[1]))
   return var_iterator

  return aggregator,annotated_var_iter

def _eval_aggregate(true_variants, pred_variants, known_fp, evaluate):
  """Evaluate by chromosome and aggregate the results."""
  dicts = []
  for chrom in true_variants.chroms.union(pred_variants.chroms):
    chrom_true = true_variants.on_chrom(chrom)
    chrom_pred = pred_variants.on_chrom(chrom)
    chrom_known_fp = known_fp.on_chrom(chrom) if known_fp else None
    dicts.append(evaluate(chrom_true, chrom_pred,chrom_known_fp))
  return _aggregate(dicts)  

def evaluate_variants(true_variants,pred_variants,eps,eps_bp,ref,window,known_fp=None):
    def evaluate(true,pred,known_fp):
      return chrom_evaluate_variants(true,pred,eps,eps_bp,ref,window,known_fp)
    return _eval_aggregate(true_variants,pred_variants,known_fp,evaluate)

# we assume that contigs are in the order provided in the .fai index
# but that contigs may be missing
def record_by_chrom(vcf_iter):
  for chrom,chrom_group in groupby(vcf_iter, lambda r: r.CHROM):
    yield chrom,chrom_group

def evaluate_low_memory(true_iter,pred_iter,eps,eps_bp,ref,window,max_indel_len,contig_lookup,writer=None,known_fp=None):
  genome_stats = chrom_var_stats_dict()
  def evaluate_low_memory_chrom(chrom_name,true_chrom,pred_chrom,known_fp=None):
    trueChromVariants = ChromVariants(chrom_name,max_indel_len)
    for r in true_chrom:
      trueChromVariants.add_record(r)
    predChromVariants = ChromVariants(chrom_name,max_indel_len)
    for r in pred_chrom:
      predChromVariants.add_record(r)
    cvs = chrom_evaluate_variants(trueChromVariants,predChromVariants,eps,eps_bp,ref,window,known_fp)
    aggregate_stats(genome_stats,cvs) # side effects, sorry
    if writer:
      write_annotated_var(writer,cvs)

  # i already hate this code.
  true_generator = record_by_chrom(true_iter)
  pred_generator = record_by_chrom(pred_iter)
  (tchrom,tchrom_records) = next(true_generator,(None,None))
  (pchrom,pchrom_records) = next(pred_generator,(None,None))

  def eval_unequal_chrom(true_tuple,pred_tuple,known_fp_tuple):
    min_chrom = [true_tuple[0],pred_tuple[0],known_fp_tuple[0]]
    def min_or_empty(t):
      if t[0] == min_chrom:
        return t[1]
      else:
        return []
    evaluate_low_memory_chrom(min_chrom,*map(min_or_empty,[true_tuple,pred_tuple,known_fp_tuple]))

  while (tchrom is not None) or (pchrom is not None):
    if tchrom == pchrom:
      evaluate_low_memory_chrom(tchrom,tchrom_records,pchrom_records)
      (tchrom,tchrom_records) = next(true_generator,(None,None))
      (pchrom,pchrom_records) = next(pred_generator,(None,None))
    elif tchrom is None:
      evaluate_low_memory_chrom(pchrom,[],pchrom_records)
      for (pchrom,pchrom_records) in pred_generator:
        evaluate_low_memory_chrom(pchrom,[],pchrom_records)
      break
    elif pchrom is None:
      evaluate_low_memory_chrom(tchrom,tchrom_records,[])
      for (tchrom,tchrom_records) in true_generator:
        evaluate_low_memory_chrom(tchrom,tchrom_records,[])
      break
    elif contig_lookup[tchrom] > contig_lookup[pchrom]:
      evaluate_low_memory_chrom(pchrom,[],pchrom_records)
      (pchrom,pchrom_records) = next(pred_generator,(None,None))
    else: # tchrom must be smaller than pchrom
      evaluate_low_memory_chrom(tchrom,tchrom_records,[])
      (tchrom,tchrom_records) = next(true_generator,(None,None))

  def get_var_stats(vartype):
    return genome_stats[vartype]
  return get_var_stats

# overwrite this. use PyVCF writer to keep all other record info
def output_annotated_variants(var_aggregate,contig_ordering,outVCF,vcf_header_lines):
    # write the header
    outVCF.write("##fileformat=VCFv4.1\n")
    for s in vcf_header_lines:
      outVCF.write(s + "\n")
    outVCF.write("##INFO=<ID=smash_type,Type=String,Description=\"classify variant as TP,FP,FN,or rescued\">\n")
    outVCF.write("##INFO=<ID=source_file,Type=Integer,Description=\"variant originally in first or second vcf passed to SMaSH\">\n")
    outVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for rec in var_aggregate(contig_ordering):
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

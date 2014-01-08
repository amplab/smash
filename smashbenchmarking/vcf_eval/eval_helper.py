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

"""Helpers for VCF evaluation."""

from __future__ import print_function

from chrom_variants import VARIANT_TYPE,GENOTYPE_TYPE,extract_range
from rectify_seq import SequenceRescuer

def _type_dict(default=0):
    return dict(map(lambda t: [t,default],VARIANT_TYPE))

def _genotype_concordance_dict():
    # this looks too much like Clojure
    gcDict = dict()
    for vtype in VARIANT_TYPE:
        gcDict[vtype] = dict()
        for gtype1 in GENOTYPE_TYPE:
            gcDict[vtype][gtype1]= dict()
            for gtype2 in GENOTYPE_TYPE:
                gcDict[vtype][gtype1][gtype2] = 0
    return gcDict

def rescue_mission(false_negatives,false_positives,loc,ref,window):
  # the rescue mission attempts to rescue a variant at a specific location in
  # the false negative track.    # note that if a variant is missed (false_negative) due to representation
  # then there will be a nearby false_positive. Thus it suffices to only
  # attempt to rescue the false negatives.
  assert false_negatives.all_variants[loc] != None # will raise a key error if broken


  # if nothing was rescued, return no new TP, no FP removed
  num_new_tp = _type_dict()
  num_fp_removed = _type_dict()
  if false_negatives.all_variants[loc].var_type.startswith("SV"):
      # don't try to rescue SVs (if they could've been rescued, they were marked true
      # by having some breakpoint within the window
      return num_new_tp,num_fp_removed
  rescuer = SequenceRescuer(false_negatives.chrom,loc,false_negatives,false_positives,ref,window)
  if not rescuer or not rescuer.rescued:
    return num_new_tp,num_fp_removed

  # now the whole truth window becomes true positives
  # and the whole predicted window is removed from false positives

  for variant in rescuer.truthWindowQueue[rescuer.windowsRescued[0]]:
    false_negatives._remove_variant(variant.pos)
    num_new_tp[variant.var_type] += 1
  for variant in rescuer.predictWindowQueue[rescuer.windowsRescued[1]]:
    false_positives._remove_variant(variant.pos)
    num_fp_removed[variant.var_type] += 1
  return num_new_tp,num_fp_removed

def var_match_at_loc(true_variants,pred_variants,loc):
    def get_var(variants):
        return variants.all_variants[loc]
    true_var = get_var(true_variants)
    pred_var = get_var(pred_variants)
    if ( true_var.var_type != pred_var.var_type ):
        return False
    return true_var.alt == pred_var.alt

def any_var_match_at_loc(fp_variants,pred_variants,loc):
    def get_var(variants):
        return variants.all_variants[loc]
    fp_var = get_var(fp_variants)
    pred_var = get_var(pred_variants)
    for a in fp_var.alt:
      for v in pred_var.alt:
        if v == a:
          return True
    return False

def allele_match_at_loc(true_variants, pred_variants, loc):
  def get_var(variants):
    return variants.all_variants[loc]
  true_var = get_var(true_variants)
  pred_var = get_var(pred_variants)
  return allele_match(true_var, pred_var)


def allele_match(true_var, pred_var):
  def alt_alleles(var):
    return set(var.alt)
  return alt_alleles(true_var) == alt_alleles(pred_var)

def structural_match(true_variant,pred_vars_all,sv_eps,sv_eps_bp):
  matches = find_possible_matches(true_variant,pred_vars_all,sv_eps_bp,sv_eps)
  if matches:
   matches = get_closest(true_variant,matches)
  return matches

class ChromVariantStats:

  """Stats for a certain contig's worth of variants."""
  
  def __init__(self, true_var, pred_var, true_positives,
               false_positives, false_negatives,concordance):
    self.chrom = true_var.chrom
    self.true_var = true_var
    self.pred_var = pred_var
    self.num_true = dict(map(lambda t: [t,true_var.var_num(t)],VARIANT_TYPE))
    self.num_pred = dict(map(lambda t: [t,pred_var.var_num(t)],VARIANT_TYPE))
    self.num_tp = _type_dict() # true positives as int
    self.num_fp = _type_dict()
    self.num_fn = _type_dict()
    self.false_positives = self._extract(pred_var,false_positives,False) #chromvariants
    self.false_negatives = self._extract(true_var,false_negatives,True) #chromvariants
    self.intersect_bad = None # set externally
    self.known_fp = None # set externally
    self.calls_at_known_fp = None # set externally
    self.known_fp_variants = None # set externally
    self.genotype_concordance = concordance
    for loc in true_positives:
      var = true_var.all_variants[loc]
      self.num_tp[var.var_type] += 1

  def _extract(self,chromvariant,locset,isFN):
    clone = chromvariant.clone()
    assert clone.all_variants == chromvariant.all_variants
    assert clone.all_locations == chromvariant.all_locations
    falseTypeDict = self.num_fn if isFN else self.num_fp if isFN != None else _type_dict()
    for loc in chromvariant.all_locations:
      if ( loc not in locset ):
        clone._remove_variant(loc)
      else:
        var = clone.all_variants[loc]
        falseTypeDict[var.var_type] += 1

    return clone

  def rectify(self, ref, window):
    """Rescue variants from VCF ambiguity.

    Given reference genome and window of sequence comparison,
    fix each error estimated to be an artifact of VCF ambiguity.

    Calling this function strictly improves realism of evaluation,
    but doesn't perfectly detect ambiguity.
    TODO: Don't do redundant calculation for overlapping windows
    TODO: This should return two bools: events match, and genotypes match

    Note: Here the window is a single integer for the size. Here it gets
    converted specifically to an interval, and that's its subsequent usage in the stack.
    """

    locs_to_rescue = list(map(lambda loc: loc, self.false_negatives.all_locations))
    # note we needed to force a copy here, since rescue_mission is modifying the false-negative sets

    for loc in locs_to_rescue:
        if ( loc in self.false_negatives.all_variants ): # if the element is still in the set of false negatives
            new_tp,rm_fp = rescue_mission(self.false_negatives,self.false_positives,loc,ref,window)
            for t in VARIANT_TYPE:
              # seemingly odd accounting. The number of predicted variants *changes* as a result of rescuing.
              # e.g. 2 predicted FPs are in fact 1 FN. So
              #  -- remove 2 predicted variants
              #  -- remove 2 false positives
              #  -- remove 1 false negative
              #  -- add 1 true positive
              self.num_pred[t] -= rm_fp[t]
              self.num_fp[t] -= rm_fp[t]
              self.num_pred[t] += new_tp[t]
              self.num_fn[t] -= new_tp[t]
              self.num_tp[t] += new_tp[t]

  def _nrd_counts(self,var_type):
    genoGenoCounts = self.genotype_concordance[var_type]
    nWrong = genoGenoCounts[GENOTYPE_TYPE.HOM_REF][GENOTYPE_TYPE.HET]
    nWrong += genoGenoCounts[GENOTYPE_TYPE.HOM_REF][GENOTYPE_TYPE.HOM_VAR]
    nWrong += genoGenoCounts[GENOTYPE_TYPE.HET][GENOTYPE_TYPE.HOM_REF]
    nWrong += genoGenoCounts[GENOTYPE_TYPE.HET][GENOTYPE_TYPE.HOM_VAR]
    nWrong += genoGenoCounts[GENOTYPE_TYPE.HOM_VAR][GENOTYPE_TYPE.HOM_REF]
    nWrong += genoGenoCounts[GENOTYPE_TYPE.HOM_VAR][GENOTYPE_TYPE.HET]
    nTotal = nWrong + genoGenoCounts[GENOTYPE_TYPE.HET][GENOTYPE_TYPE.HET] + genoGenoCounts[GENOTYPE_TYPE.HOM_VAR][GENOTYPE_TYPE.HOM_VAR]
    return (nWrong,nTotal)

  def to_dict(self,var_type):
    stats = {}
    stats['num_true'] = self.num_true[var_type]
    stats['num_pred'] = self.num_pred[var_type]
    stats['false_positives'] = self.num_fp[var_type]
    stats['false_negatives'] = self.num_fn[var_type]
    stats['good_predictions'] = self.num_tp[var_type]
    stats['intersect_bad'] = len(self.intersect_bad[var_type])
    nrd_wrong,nrd_total = self._nrd_counts(var_type)
    stats['nrd_wrong'] = nrd_wrong
    stats['nrd_total'] = nrd_total
    stats['known_fp_calls'] = self.calls_at_known_fp[var_type] if self.calls_at_known_fp else 0
    stats['known_fp'] = self.known_fp[var_type] if self.calls_at_known_fp else 0
    return stats

def chrom_evaluate_variants(true_var,pred_var,known_fp,sv_eps,sv_eps_bp,ref,window):
    true_loc = set(true_var.all_locations)
    pred_loc = set(pred_var.all_locations)
    genotype_concordance = _genotype_concordance_dict()
    if pred_var.negative_snp_labels:
        raise Exception("Why do you have monomorphic calls?")
    neg = true_var.negative_snp_labels
    false_positives = pred_loc.difference(true_loc)
    false_negatives = true_loc.difference(pred_loc)
    intersect_good = []
    intersect_bad = []
    intersect_bad_dict = _type_dict([])
    for loc in pred_loc.intersection(true_loc):
      vartype = true_var.all_variants[loc].var_type
      if ( vartype.startswith("SV") ): # ignore SVs here
        continue
      match = var_match_at_loc(true_var, pred_var, loc)
      destination = intersect_good if match else intersect_bad
      destination.append(loc)
      if not match:
       intersect_bad_dict[vartype].append(loc)
      else:
       true_geno = true_var.all_variants[loc].genotype_type
       pred_geno = pred_var.all_variants[loc].genotype_type
       genotype_concordance[vartype][true_geno][pred_geno] += 1

    # match also calls at known false positives
    known_fp_calls_positions = list()
    calls_at_known_fp = _type_dict()
    all_known_fp = _type_dict()
    if known_fp:
      for loc in pred_loc.intersection(known_fp.all_locations):
        vartype = known_fp.all_variants[loc].var_type
        match = any_var_match_at_loc(known_fp,pred_var,loc)
        if match:
         calls_at_known_fp[vartype] += 1
         known_fp_calls_positions.append(loc)
        all_known_fp[vartype] += 1

    # structural variants are a special case if not matching exactly
    for loc in (true_loc - pred_loc):
      vartype = true_var.all_variants[loc].var_type
      if ( not vartype.startswith("SV") ):
        continue
      match = structural_match(true_var.all_variants[loc],pred_var,sv_eps,sv_eps_bp)
      if match and match in false_positives:   # don't double count
        intersect_good.append(loc)
        false_positives.remove(match)
        false_negatives.remove(loc)
        true_geno = true_var.all_variants[loc].genotype_type
        pred_geno = pred_var.all_variants[match].genotype_type
        genotype_concordance[vartype][true_geno][pred_geno] += 1

    false_positives.update(intersect_bad)
    false_negatives.update(intersect_bad)
    variant_stats = ChromVariantStats(true_var, pred_var,
                                      intersect_good, false_positives,
                                      false_negatives,genotype_concordance)
    if ( window and ref):
      variant_stats.rectify(ref, window)

    if ( known_fp ):
      variant_stats.known_fp = all_known_fp
      variant_stats.calls_at_known_fp = calls_at_known_fp
      variant_stats.known_fp_variants = variant_stats._extract(variant_stats.pred_var,known_fp_calls_positions,None)

    variant_stats.intersect_bad = intersect_bad_dict
    #stats = variant_stats.to_dict()
    #stats['intersect_bad'] = len(intersect_bad)
    return variant_stats

def is_within(low, mid, high, eps):
  """General helper."""
  return low - eps <= mid <= high + eps


def indel_or_sv_match(true_var, pred_var, eps_bp, eps_len):
  """Does a true INDEL/SV match a predicted one, at given tolerances?

  Yes iff:
    a) predicted and true variant have the same type,
       e.g. insertion or deletion
    b) predicted breakpoint within 'eps_bp' of true variant,
       allowing for ambiguity due to "bookending"
    c) length of deletion/insertion within 'eps_len' of true length

  Don't check zygosity; for heterozygous variants, allow any choice of allese.
  """
  # TODO (post paper) benchmark zygosity of indels and possibly SVs
  # not all callers report it
  # TODO if ignoring genotype/benchmarking variant discovery, check that set of gains is the same
  # TODO (post paper) check sequence for indels
  true_type = true_var.var_type
  pred_type = pred_var.var_type
  assert 'SNP' not in true_type and 'SNP' not in pred_type
  true_pos = true_var.pos

  pos_ok = is_within(true_pos, pred_var.pos, true_pos, eps_bp)
  def gains_ok(true, pred):
    return is_within(true, pred, true, eps_len)
  some_gains_ok = any(gains_ok(true, pred) \
                        for true in true_var.gains for pred in pred_var.gains)
  return true_type == pred_type and pos_ok and some_gains_ok


def find_possible_matches(var, pred_variants, eps_bp, eps_len):
  """Find all predicted locations that could possibly match a true location.

  @param var: a true variant
  """
  min_loc = var.pos - eps_bp
  max_loc = var.pos + eps_bp
  var_type = var.var_type
  candidate_locs = extract_range(pred_variants.var_locations(var_type),
                                 min_loc, max_loc)
  def is_match(loc):
    pred_var = pred_variants.var(var_type, loc)
    return indel_or_sv_match(var, pred_var, eps_bp, eps_len)
  return filter(is_match, candidate_locs)


def get_closest(true_var, possible_matches):
  """Return the possible match closest to a given variant.

  If there are multiple closest, pick the first one.
  """
  def get_distance(pred_loc):   # TODO: Maybe use length or junction.
    return abs(true_var.pos - pred_loc)
  min_distance = min(map(get_distance, possible_matches))
  def is_closest(possible_match):
    return get_distance(possible_match) == min_distance
  closests = filter(is_closest, possible_matches)
  return closests[0]
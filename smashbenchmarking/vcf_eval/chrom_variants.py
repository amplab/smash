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

"""Variation on a single chromosome."""

from __future__ import print_function
import vcf
from collections import defaultdict
import bisect
import sys

class Enum(set):
    def __getattr__(self,name):
        if name in self:
            return name
        raise AttributeError(name)

VARIANT_TYPE = Enum(["SNP","INDEL_DEL","INDEL_INS","INDEL_OTH","SV_DEL","SV_INS","SV_OTH"])
GENOTYPE_TYPE = Enum(["HOM_REF","HET","HOM_VAR","NO_CALL","UNAVAILABLE"])

# this maps the PyVCF representation to this enum representation
GENOTYPE_TYPE_MAP = { 0 : GENOTYPE_TYPE.HOM_REF, 1 : GENOTYPE_TYPE.HET, 2 : GENOTYPE_TYPE.HOM_VAR,
                      None:GENOTYPE_TYPE.NO_CALL}

def _lacks_alt_alleles(record):
  alt = map(str, record.ALT)
  return alt[0] == 'None' or not alt[0]


def is_pass(record):
  """Return whether a PyVCF record passes all filters.

  From the VCF spec, that's if it's 'PASS' or '.'.
  But PyVCF uses '[]' for 'PASS' (older versions used 'None').
  PyVCF uses None for "." now, presumably to be pythonic (not [] is true; not None is true)
  """
  filters = record.FILTER
  return not filters or filters == '.'

def is_sv(record,maxSize):
    """
    return whether a PyVCF record is an SV or not (by our definition). This is either of
    size > maxSize
    symbolic alleles
    @param record:  the record
    @return:
    """
    refSize = len(record.REF)
    altSizes = map(len,filter(lambda a: not isinstance(a, vcf.model._SV) and a is not None, record.ALT))
    symbolic = any(map(lambda t: "<" in t,str(record.ALT)))
    if ( refSize > maxSize ):
        return True
    if ( any(map(lambda t: t>maxSize,altSizes)) ):
        return True
    return symbolic

class ChromVariants:

  """Variation on a single chromosome."""

  def __init__(self, chrom, max_indel_len, **kwargs):
    """All variants in a VCF file, parsed and classified.
    """
    self._args = kwargs
    self.chrom = chrom
    # Hold variants by type; list vs dict
    self._var_locations = defaultdict(list)
    self.snp_pos_dict = {} 
    self.indel_pos_dict = {} 
    self.sv_pos_dict = {} 
    self._max_indel_len = max_indel_len
    # Hold variants not segregated by type; list vs dict
    self.all_locations = []
    self.all_variants = {}
    self.negative_snp_labels = set() # For Sampled Human datasets of SMaSH.

  def clone(self):
    newVar = ChromVariants(self.chrom,self._max_indel_len,**self._args)
    for loc in self.all_locations:
        v = self.all_variants[loc]
        newVar._add_variant(v)
    return newVar

  def _var_iterator(self,locs=None):
      if ( locs == None ):
          locs = self._var_locations
      return map(lambda t: self.all_variants[t],locs)

  def _vcf_iterator(self,info=".",locs=None):
      return map(lambda s: s._vcf_entry(self.chrom,info), self._var_iterator(locs)).__iter__()

  def _add_variant(self,var):
    if ( var.pos in self.all_variants ):
        #raise Exception("Two variants have the exact same start position: "+ str(self.all_variants[var.pos]) + 'and' + str(var))
        print("Two variants have the exact same start position: "+ str(self.all_variants[var.pos]) + ' and ' + str(var)+". Variant excluded.",file=sys.stderr )
        return
    # careful as adding a variant requires a later sorting of all the positions
    self._var_locations[var.var_type].append(var.pos)
    self._var_dict(var.var_type)[var.pos] = var
    self.all_locations.append(var.pos)
    self.all_variants[var.pos] = var

  def _ensure_sorted(self):
      for t in VARIANT_TYPE:
          self._var_locations[t].sort()
      self.all_locations.sort()

  def _remove_variant(self,loc):
      var = self.all_variants[loc]
      self._var_locations[var.var_type].remove(loc)
      self._var_dict(var.var_type).pop(loc)
      self.all_locations.remove(loc)
      self.all_variants.pop(loc)
      return var

  def validate(self):
      diff =  len(self.all_locations) - len(self.all_variants)
      if ( diff != 0 ):
          # this is bad.
          except_text = "Exception:\nThe number of variant positions does not match the number of variants.\n"
          # get some positions that aren't in the dict
          no_variants = set(self.all_locations).difference(self.all_variants.keys())
          no_positions = map(lambda t: self.all_variants[t], set(self.all_variants.keys()).difference(self.all_locations))
          except_text += "Missing variants at: \n"
          except_text += "\n".join(map(lambda i: str(no_variants[i] if i < len(no_positions) else "None"),range(3)))
          except_text += "\n Missing positions for: \n"
          except_text +=  "\n".join(map(lambda i: str(no_positions[i] if i < len(no_positions) else "None"),range(3)))
          except_text += "\n Potential duplicate elements of position list: \n"
          positions = set(self.all_locations)
          position_list = list(map(lambda a: a, self.all_locations))
          for p in self.all_locations:
              if ( p in positions ):
                  positions.remove(p)
                  position_list.remove(p)
          except_text += "\n".join(map(lambda i: str(position_list[i] if i < len(position_list) else "None"),range(3)))
          raise AssertionError(except_text)

  def _var_dict(self, var_type):
    for kind in ['SNP', 'INDEL', 'SV']:
      if kind in var_type:
        return getattr(self, kind.lower() + '_pos_dict')
    AttributeError('Bad variant type: ' + var_type)

  def add_record(self, record):

    """Parse a VCF record.

    Assumes only one 'indel_or_sv' at each location,
    but doesn't check if we also have a variant of the other two kinds.
    """

    assert record.CHROM == self.chrom
    if not is_pass(record) and not self._args.get('knownFP',None): # if this is a known FP object, let it past the filter
      return
    ref = record.REF
    if ( record.REF != record.REF.upper() ):
      raise AssertionError("VCF contains lower-case bases in the reference: " + record.REF+ " at "+str(record.POS))

    alt = map(str, record.ALT)
    if _lacks_alt_alleles(record):
      raise Exception("Monomorphic records (no alt allele) are not supported.")

    len_ref = len(ref)
    assert len(record.samples) == 1 # We expect vcfs to correspond to a
                                    # single sample.
    # Possible TODO: filter on desired sample? Or at least document that this is our convention.
    # given that there's a single sample, get the genotype type
    genotype_type = GENOTYPE_TYPE_MAP[record.samples[0].gt_type]

    if (not self._args.get('knownFP',None)) and (not is_sv(record, self._max_indel_len)) and ( genotype_type == GENOTYPE_TYPE.HOM_REF or genotype_type == GENOTYPE_TYPE.NO_CALL ):
        raise Exception("Monomorphic records (0/0 or ./. or . genotypes) are not supported (except in known FP VCFs). Observed: "+str(record))

    def add_variant(var_type):
      pos = record.POS
      new_variant = Variant(pos, ref, alt, var_type, genotype_type)
      self._add_variant(new_variant)
      

    def add_indel_or_sv(is_indel):
      def add_appropriate_variant(indel_type,sv_type):
        if ( is_indel ):
            add_variant(indel_type)
        else:
            add_variant(sv_type)

      if len(alt) == 1 and len(alt[0]) == 1:
        add_appropriate_variant(VARIANT_TYPE.INDEL_DEL,VARIANT_TYPE.SV_DEL)
      elif len_ref == 1:
        add_appropriate_variant(VARIANT_TYPE.INDEL_INS,VARIANT_TYPE.SV_INS)
      else:
        add_appropriate_variant(VARIANT_TYPE.INDEL_OTH,VARIANT_TYPE.SV_OTH)
        # Other now contains MNP, inversions
        # MNP question: can ref/alt be of different lengths?
        # MNP: all indels that are not currently being handled?
        # Short inversions may be MNPs
        # Inversion: ACGT -> ATGC
        # This is all dependent on our validation data.

        # MNP: many to many shorter than indel constant
        # inversion: many to many longer than indel constant that are reversed
        # all others: still falls in other buckets

    if record.is_snp:
      add_variant(VARIANT_TYPE.SNP)
    else:
      allele_lengths = map(len, alt) + [len_ref]
      is_indel = (not is_sv(record,self._max_indel_len) ) and max(allele_lengths) <= self._max_indel_len
      add_indel_or_sv(is_indel)

  def var_locations(self, var_type):
    return self._var_locations[var_type]

  def var_num(self, var_type):
    return len(self.var_locations(var_type))

  def var(self, var_type, loc):
    """Lookup a variant given type and location."""
    return self._var_dict(var_type)[loc]

  def getAllVariants(self):
      return map(lambda t: self.all_variants[t],self.all_locations)


class Variant:

  def __init__(self, pos, ref, alt, var_type, genotype_type):
    self.__mutable = True
    self.pos = pos 
    self.ref = ref 
    self.alt = alt 
    self.var_type = var_type
    self.genotype_type = genotype_type

  @property
  def gains(self):
      return map(lambda allele: len(allele) - len(self.ref), self.alt)

  @property
  def losses(self):
      return map(lambda g: max(0,-g),self.gains)

  def __str__(self):
      return "%d %s/%s %s %s" % (self.pos,self.ref,str(self.alt),self.var_type,self.genotype_type)

  def overlaps_allele(self,pos):
      return any(map(
          lambda loss: self.pos + loss >= pos and self.pos <= pos, # can only overlap in the forward direction
          self.losses))

  def overlaps(self,pos):
      return self.pos + len(self.ref) >= pos and self.pos <= pos

  def strictly_overlaps(self,pos):
      return self.pos + len(self.ref) > pos and self.pos <= pos

  def _vcf_entry(self,contig,info="."):
      return "%s\t%d\t.\t%s\t%s\t.\t.\t%s" % (contig,self.pos,self.ref,",".join(self.alt),info)

def _getRestOfPath(chosenSoFar,remainingChoices):
    if ( remainingChoices == [] ):
        return [chosenSoFar]
    choices = filter(lambda choice: not any(map(lambda chosen: chosen.strictly_overlaps(choice.pos),chosenSoFar)),remainingChoices[0])

    def getRest(c):
        nowChosen = chosenSoFar + [c]
        return _getRestOfPath(nowChosen,remainingChoices[1:len(remainingChoices)])

    # this recursion returns nested paths down the decision tree
    # the control function flattens subpaths where there are choices
    control = lambda retpaths: [choice for path in retpaths for choice in path] # flatten

    return control(map(getRest,choices))

def _getOverlaps(overlapSets,variants):
    if ( variants == [] ):
            return overlapSets

    if ( overlapSets == [] ):
        overlapSets.append([variants[0]])
        return _getOverlaps(overlapSets,variants[1:len(variants)])

    nextVar = variants[0]
    if ( any(map(lambda var: var.strictly_overlaps(nextVar.pos),overlapSets[len(overlapSets)-1])) ):
        overlapSets[len(overlapSets)-1].append(nextVar)
    else:
        overlapSets.append([nextVar])
    return _getOverlaps(overlapSets,variants[1:len(variants)])

def extract_range(pred_locs, min_loc, max_loc):
  """General helper: restrict a sorted list to a range.
   @pred_locs - a list of positions SORTED IN ASCENDING ORDER
  """
  low = bisect.bisect_left(pred_locs, min_loc)
  high = bisect.bisect_right(pred_locs, max_loc-1)
  return pred_locs[low:high]

def extract_range_and_filter(variants,min_loc,max_loc,loc_of_interest):
  """ Given a range of locations (min_loc, max_loc), a list of variant locations (variant_locs),
      and a collection of variants (variants.all_variants), restrict the collection of variants only
      to those falling within (min_loc,max_loc).
      Then, given a loc_of_interest, remove any variant that overlaps the loc of interest.
  """
  locs_in_window = extract_range(variants.all_locations,min_loc,max_loc)
  variants_in_window = map(lambda t: variants.all_variants[t],locs_in_window)
  variants_in_window = filter(lambda var: not var.var_type.startswith("SV"),variants_in_window) # don't include SVs
  if ( loc_of_interest not in locs_in_window ):
    # trying to rescue something in the other track, no need to do any removal
    return variants_in_window
  # now remove anything that overlaps the variant at the loc of interest
  voi_size = max(variants.all_variants[loc_of_interest].losses)
  variants_in_window = filter(lambda var:
                        var.pos == loc_of_interest or # keep the loc of interest or anything not overlapping with it
                        not (var.overlaps_allele(loc_of_interest) or var.overlaps_allele(loc_of_interest+voi_size)),variants_in_window )
  return variants_in_window


# def verifyPaths(variantPaths):
#     # verify that no overlapping variants are included. This should be removed from the stack after testing.
#     for path in variantPaths:
#         soFar = []
#         for elem in path:
#             for otherElem in soFar:
#                 if ( elem == otherElem or otherElem.overlaps_allele(elem) ):
#                     raise AssertionError("Invalid Path: "+str(map(str,path)))
#     return variantPaths

def extract_variant_queues(variants,min_loc,max_loc,loc_of_interest):
  """ Given a range of locations (min_loc, max_loc), a list of variant locations (variant_locs),
      and a collection of variants (variants.all_variants), restrict the collection of variants only
      to those falling within (min_loc,max_loc).

      Then, given a loc_of_interest, remove any variant that overlaps the loc of interest.

      Then, convert the resulting list into a collection of variants.

      If any further variants overlap, choose every possible collection of variants resulting from choosing
      one and only one of the overlapping variants.
  """
  variants_in_window = extract_range_and_filter(variants,min_loc,max_loc,loc_of_interest)
  if not variants_in_window:
    return []

  # get the overlap sets TODO: inefficient - variant size can be used intelligently
  overlaps_sets = _getOverlaps([],variants_in_window)

  # return all the possible paths through the overlap sets
  return _getRestOfPath([],overlaps_sets) # verifyPaths(_getRestOfPath([],overlaps_sets))




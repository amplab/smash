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

"""Rectify scoring errors due to VCF ambiguity."""

from __future__ import print_function, division

import sys
import bisect
import chrom_variants
from chrom_variants import GENOTYPE_TYPE,VARIANT_TYPE

WINDOW_SIZE_LIMIT = 5000 # if the window size gets too big, kills not only performance but also window will likely
                         # fail to rescue
WINDOW_MAX_OVERLAPPING = 16 # if there are enough overlapping variants in the window to generate
                          # this many (or more) sequences, don't try to attempt to rescue the variants

WINDOW_VARIANT_LOOKBACK_SIZE = 50 # this should match max_indel_length

class ReferenceMismatchError(Exception):
    def __init__(self, msg):
        self.msg = msg

class SequenceRescuer(object):

    def __init__(self,contig,location,falseNegatives,falsePositives,truePositives,reference,windowSize):
        """
        @params contig: String name of chromosome
        @params location: int location of variant to rescue
        @params falseNegatives: ChromVariants holding false negatives
        @params falsePositives: ChromVariants holding false falsePositives
        @params reference: Genome for reference
        @params windowSize: int size of window for rescue
        """
        self.contig = contig
        self.window = self._getWindowSize(location,falseNegatives,falsePositives,truePositives,windowSize)
        if ( self.window[1]-self.window[0] > WINDOW_SIZE_LIMIT ):
            # window is too big; abort
            self.rescued = False
            return

        self.truthWindowQueue = chrom_variants.extract_variant_queues(falseNegatives,self.window[0],self.window[1]-1,location)
        self.predictWindowQueue = chrom_variants.extract_variant_queues(falsePositives,self.window[0],self.window[1]-1,location)

        self.truePositives = chrom_variants.extract_range_and_filter(truePositives,self.window[0],self.window[1]-1,location)

        if ( not self.truthWindowQueue or not self.predictWindowQueue or len(self.truthWindowQueue) * len(self.predictWindowQueue) > WINDOW_MAX_OVERLAPPING ):
            # either no variants or too many overlapping variants to check; abort
            self.rescued = False
            return
        try:
            self.rescued, self.windowsRescued = self._try_rescue(reference)
            self.rescued_GA = False # self.rescued and self._try_rescue_window(reference,self.windowsRescued[0],self.windowsRescued[1],True)
            #NB: rescued_GA property is never used
        except ReferenceMismatchError as err:
            print(err.msg,file=sys.stderr)
            self.rescued = False
            return
        except AssertionError:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            import traceback
            errorText = "Error during rescue attempt. Attempting to rescue variant "
            errorText += str(falseNegatives.all_variants[location])
            errorText += " at location "+ str(location)
            errorText += " with window "+str(self.window)
            errorText += "\n. The truth window queue was: \n\n"
            errorText += str(map(lambda a: map(str,a),self.truthWindowQueue))
            errorText += "\n\n And the predicted window queue was: \n\n"
            errorText += str(map(lambda a: map(str,a),self.predictWindowQueue))
            errorText += "\n"
            errorText += "This can happen if the lookback for variants overlapping a window is too small, and you have many large deletions."
            errorText += "\nConsider manually increasing the lookback in get_chopped_variants within this file.\n\n Caught assertion: \n\n"
            errorText += "".join(traceback.format_exception(exc_type,exc_value,exc_traceback))
            raise AssertionError(errorText)

    def _getWindowSize(self,location,fnVar,fpVar,tpVar,size):
        low,high = location-size,location+size
        prevL,prevH = None,None
        while ( low != prevL or high != prevH ):
            prevL,prevH = low,high
            low,high = _enlarge_bounds(fnVar,low,high)
            low,high = _enlarge_bounds(fpVar,low,high)
            low,high = _enlarge_bounds(tpVar,low,high)
            assert low <= prevL and high >= prevH
        return (low,high,self.contig)

    def _try_rescue(self,ref):
        for true_idx in range(len(self.truthWindowQueue)):
            for pred_idx in range(len(self.predictWindowQueue)):
                # we always ignore genotype for now
                if ( self._try_rescue_window(ref,true_idx,pred_idx,False) ):
                    return True, (true_idx,pred_idx)
        return False,None

    def _add_true_pos_to_queue(self,queue):
        new_queue = list(queue) # make copy as we're altering in place
        for tp in self.truePositives:
            for var in queue:
                if tp.strictly_overlaps_var(var):
                    return False # if a true pos overlaps with a variant in the queue, don't rescue
            new_queue.append(tp)
        new_queue.sort(key=lambda v: v.pos)
        return new_queue

    def _try_rescue_window(self,ref,window_true_num,window_pred_num,genotypeAware):
        true_vars = self.truthWindowQueue[window_true_num]
        pred_vars = self.predictWindowQueue[window_pred_num]

        # at least one must have something which is not a SNP in order to rescue
        has_non_snp = reduce(lambda a,b: a or b,map(lambda u: u.var_type != VARIANT_TYPE.SNP, true_vars+pred_vars))

        if ( not has_non_snp ):
            return False

        fn_tp_vars = self._add_true_pos_to_queue(true_vars)
        fp_tp_vars = self._add_true_pos_to_queue(pred_vars)

        if (not fn_tp_vars or not fp_tp_vars):
            print("failed when checking new queues")
            print("fn tp vars " + str(fn_tp_vars))
            return False # true pos overlapped with something in queue so we can't rescue

        trueSeq,trueSeqHet = _get_seq(self.window,fn_tp_vars,ref,genotypeAware)
        predSeq,predSeqHet = _get_seq(self.window,fp_tp_vars,ref,genotypeAware)

        return (trueSeq == predSeq) and ( trueSeqHet == predSeqHet ) # second statement always true if not genotype aware




# General helpers.

def print_debug(message):
  # print("DEBUG:", message, file=sys.stderr)
  pass

def find_le(a, x):
  """Find rightmost value less than or equal to x.

  General helper, adapted from <http://docs.python.org/2/library/bisect.html#searching-sorted-lists>.
  """
  i = bisect.bisect_right(a, x)
  if i:
    return a[i - 1]
  raise ValueError


def _get_seq(window,variants,ref,genotypeAware):
    """
    Using the variation in @variants, construct two haplotypes, one which
    contains only homozygous variants, the other which contains both hom and het variants
    by placing those variants into the reference base string
    @param variants: A vcf_eval.ChromVariants object
    @param low: the starting position
    @param high: the ending position
    @param ref: a parsers.genome object
    @param loc: the location that we are trying to rescue
    @param genotype: whether to phase hets onto their own sequence to check for genotype accuracy (if there are multiple and they don't overlap, phasing doesn't matter)
    @return: a tuple of sequences of bases that comes from modifying the reference sequence with the variants
    """
    low = window[0]
    high = window[1]
    hetChunks = []
    homChunks = []
    hetOffset = low
    homOffset = low

    # note: if genotypeAware is False, the het chunks/offset will not be used

    def get_ref_bases(start,end):
        """VCF parser is 1-based, but genome is 0-based."""
        return ref.ref(window[2],start-1,end-1)
    def add_ref_bases_until(chunks,begin,end):
        chunks.append(get_ref_bases(begin,end))
    def add_alt(chunk,start,var):
        add_ref_bases_until(chunk,start,var.pos)
        chunk.append(var.alt[0])

    for variant in variants:
        loc = variant.pos
        #print((variant.ref, get_ref_bases(variant.pos,variant.pos+len(variant.ref))))
        verifyRefBases = get_ref_bases(variant.pos,variant.pos+len(variant.ref))
        if ( variant.ref != verifyRefBases ):
            raise ReferenceMismatchError("Variant ref does not match reference at " + window[2] + " " + str(loc) + ": " +variant.ref + " != " + verifyRefBases )
        assert hetOffset <= loc and homOffset <= loc
        assert variant.genotype_type != GENOTYPE_TYPE.HOM_REF
        assert variant.genotype_type != GENOTYPE_TYPE.NO_CALL
        if ( (not genotypeAware) or variant.genotype_type == GENOTYPE_TYPE.HOM_VAR):
            add_alt(homChunks,homOffset,variant)
            homOffset = len(variant.ref) + loc
        else: # ( variant.genotype_type == GENOTYPE_TYPE.HET )
            add_alt(hetChunks,hetOffset,variant)
            hetOffset = len(variant.ref) + loc
# NB: this check seems redundant with the assert after it
        if ( hetOffset > high or homOffset > high ):
            print("-----fail-----")
            print(window)
            print(map(str,variants))
            print((homOffset,high))
        assert hetOffset <= high and homOffset <= high
    if ( genotypeAware ):
        add_ref_bases_until(hetChunks,hetOffset,high)
    add_ref_bases_until(homChunks,homOffset,high)
    return (''.join(homChunks),''.join(hetChunks))

def _get_chopped_variant(variants, loc, higher):
  """If 'loc' points to a base in a ref allele, return the variant it chops.

  Else return 'None'. Note: now that variants are allowed to overlap it's possible to have

  V1:    XXXXXXXXXXX-----XX|XX
  V2:    XXXXXXXXX-----------X

  where V2 overlaps the window but V1 does not (though V1 starts after V2).
  """
  locations = variants.all_locations
  def getNear(p):
      if not p:
          return None
      try:
          return find_le(locations,p)
      except ValueError:
          return None

  vars_near_bound = [getNear(loc)]

  if vars_near_bound[0] is None:
    return None

  dist = loc-vars_near_bound[0]

  offset = 0
  while dist <= WINDOW_VARIANT_LOOKBACK_SIZE : # up this if we expect variants further away to exceed the boundary (?)
      if vars_near_bound[offset]:
          dist = loc-vars_near_bound[offset]
          vars_near_bound.append(getNear(vars_near_bound[offset]-1))
          offset += 1
      else:
          break

  # remove the Nones
  vars_near_bound = map(lambda t: variants.all_variants[t],filter(lambda u: u != None, vars_near_bound))

  if ( not vars_near_bound ):
      return None

  assert not any(map(lambda t: t.pos > loc,vars_near_bound))

  overlapping_vars_near_bound = filter(lambda t: t.overlaps(loc), vars_near_bound)
  if ( not overlapping_vars_near_bound ):
    return None

  if higher:
    overlapping_vars_near_bound.sort(key=lambda t: (t.pos + len(t.ref)), reverse=True)
  else:
    overlapping_vars_near_bound.sort(key=lambda t: t.pos)

  return overlapping_vars_near_bound[0]



def _enlarge_bounds(variants, low, high):
  """Enlarge the range [low, high) to not chop variants.

  Assume variants aren't next to each other.
  """
  variant = _get_chopped_variant(variants, low, False)
  if variant:
    low = variant.pos
  variant = _get_chopped_variant(variants, high - 1,True)
  if variant:
    if ( high == variant.pos + len(variant.ref) ):
        # edge case - variant abuts the window
        high += 1
    elif ( high < variant.pos + len(variant.ref) ):
        high = variant.pos + len(variant.ref)

  return low, high

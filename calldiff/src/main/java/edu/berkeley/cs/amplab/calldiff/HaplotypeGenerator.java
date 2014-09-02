package edu.berkeley.cs.amplab.calldiff;

import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The algorithm for generating all the possible haplotypes over a genomic window. The calls are
 * assumed to all be non-overlapping. First, the calls are partitioned into phasesets
 * (unphased calls are singleton phasesets). Then, for each allele, we iterate down the haplotype
 * and insert the appropriate replacement. One complicating factor is that since all call
 * coordinates are relative to the original reference sequence, for each haplotype generated we
 * need to maintain a mapping between the genomic coordinate space and the haplotype string
 * coordinate space.
 */
public class HaplotypeGenerator {

  static class Haplotype {

    static Haplotype create(String string, int genomeCoordinate) {
      NavigableMap<Integer, Integer> offsets = new TreeMap<>();
      offsets.put(genomeCoordinate, 0);
      return new Haplotype(string, offsets);
    }

    private final String string;
    private final NavigableMap<Integer, Integer> offsets;

    private Haplotype(String string, NavigableMap<Integer, Integer> offsets) {
      this.string = string;
      this.offsets = offsets;
    }

    private static <X> X getProperty(List<Call> calls, Function<? super Call, ? extends X> accessor,
        String message) {
      List<X> list = calls.stream().map(accessor).distinct().collect(Collectors.toList());
      Preconditions.checkState(1 == list.size(), message);
      return Iterables.getOnlyElement(list);
    }

    Stream<Haplotype> apply(List<Call> calls) {
      getProperty(calls, Call::phaseset, "All calls must belong to the same phaseset").orElseGet(
          () -> {
            Preconditions.checkState(1 == calls.size(), "Only 1 unphased call allowed");
            return null;
          });
      Stream.Builder<Haplotype> haplotypes = Stream.builder();
      for (int numAlleles = getProperty(calls, call -> call.genotype().size(),
          "Calls in same phaseset has different number of alleles"), allele = 0;
          allele < numAlleles; ++allele) {
        StringBuilder newString = new StringBuilder(string);
        NavigableMap<Integer, Integer> newOffsets = new TreeMap<>();
        newOffsets.putAll(offsets);
        for (Call call : calls) {
          String
              reference = call.reference(),
              replacement = Stream
                  .concat(Stream.of(reference), call.alternates().stream())
                  .collect(Collectors.toList())
                  .get(call.genotype().get(allele));
          int referenceLength = reference.length(),
              replacementLength = replacement.length(),
              delta = replacementLength - referenceLength,
              start = call.position(),
              end = start + referenceLength;
          Map.Entry<Integer, Integer> floor = newOffsets.floorEntry(start);
          int newStart = start + floor.getValue() - floor.getKey(),
              newEnd = newStart + referenceLength;
          assert newString.substring(newStart, newEnd).equals(reference);
          newString.replace(newStart, newEnd, replacement);
          NavigableMap<Integer, Integer> offsetsCopy = new TreeMap<>();
          newOffsets.subMap(Integer.MIN_VALUE, true, end, false)
              .entrySet()
              .stream()
              .forEach(entry -> offsetsCopy.put(entry.getKey(), entry.getValue()));
          offsetsCopy.put(end, newStart + replacementLength);
          newOffsets.subMap(end, false, Integer.MAX_VALUE, true)
              .entrySet()
              .stream()
              .forEach(entry -> offsetsCopy.put(entry.getKey(), entry.getValue() + delta));
          newOffsets = offsetsCopy;
        }
        haplotypes.add(new Haplotype(newString.toString(), newOffsets));
      }
      return haplotypes.build();
    }

    @Override
    public String toString() {
      return string;
    }
  }

  public static Set<String> generateHaplotypes(FastaReader.FastaFile reference, String contig,
      List<Call> calls, int beginning, int end) {
    Stream<Haplotype> haplotypes =
        Stream.of(Haplotype.create(reference.get(contig, beginning, end), beginning));
    for (List<Call> phaseset : partitionByPhaseset(calls)) {
      haplotypes = haplotypes.flatMap(haplotype -> haplotype.apply(phaseset));
    }
    return haplotypes.map(Object::toString)
        .collect(Collectors.toSet());
  }

  static List<List<Call>> partitionByPhaseset(List<Call> calls) {
    Stream.Builder<Stream.Builder<Call>> partition = Stream.builder();
    Map<Call.Phaseset, Stream.Builder<Call>> bucketCache = new HashMap<>();
    calls.stream().forEach(call -> call.phaseset()
        .map(phaseset -> Optional.ofNullable(bucketCache.get(phaseset)).orElseGet(() -> {
              Stream.Builder<Call> bucket = Stream.builder();
              partition.add(bucket);
              bucketCache.put(phaseset, bucket);
              return bucket;
            }))
        .orElseGet(() -> {
              Stream.Builder<Call> bucket = Stream.builder();
              partition.add(bucket);
              return bucket;
            })
        .add(call));
    return partition.build()
        .map(stream -> stream.build().collect(Collectors.toList()))
        .collect(Collectors.toList());
  }
}

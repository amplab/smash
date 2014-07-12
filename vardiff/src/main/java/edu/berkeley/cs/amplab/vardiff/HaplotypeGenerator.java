package edu.berkeley.cs.amplab.vardiff;

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

public class HaplotypeGenerator {

  static class Haplotype {

    private static final Function<Call, Integer> NUMBER_OF_ALLELES = call -> call.genotype().size();
    private static final Function<Call, Optional<Call.Phaseset>> PHASESET = Call::phaseset;

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

    private static <X> X
        getProperty(List<Call> calls, Function<? super Call, ? extends X> accessor, String message) {
      List<X> list = calls.stream().map(accessor).distinct().collect(Collectors.toList());
      Preconditions.checkState(1 == list.size(), message);
      return Iterables.getOnlyElement(list);
    }

    Stream<Haplotype> apply(List<Call> calls) {
      getProperty(calls, PHASESET, "All calls must belong to the same phaseset").orElseGet(() -> {
        Preconditions.checkState(1 == calls.size(), "Only 1 unphased call allowed");
        return null;
      });
      Stream.Builder<Haplotype> haplotypes = Stream.builder();
      for (
          int numAlleles = getProperty(calls, NUMBER_OF_ALLELES, "Different numbers of alleles on calls"), allele = 0;
          allele < numAlleles;
          ++allele) {
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
          newOffsets.put(end, newStart + replacementLength);
          for (
              Map.Entry<Integer, Integer> entry = newOffsets.ceilingEntry(end + 1);
              null != entry;
              entry = newOffsets.higherEntry(entry.getKey())) {
            entry.setValue(entry.getValue() + delta);
          }
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

  public static Set<String> generateHaplotypes(
      FastaReader.FastaFile reference,
      List<Call> calls,
      int beginning,
      int end) {
    List<List<Call>> partition = partitionByPhaseset(calls);
    throw new UnsupportedOperationException();
  }

  static List<List<Call>> partitionByPhaseset(List<Call> calls) {
    Stream.Builder<Stream.Builder<Call>> partition = Stream.builder();
    Map<Call.Phaseset, Stream.Builder<Call>> bucketCache = new HashMap<>();
    calls.stream().forEach(call -> {
      call.phaseset()
          .map(phaseset -> Optional.ofNullable(bucketCache.get(phaseset))
              .orElseGet(() -> {
                Stream.Builder<Call> bucket = Stream.builder();
                bucketCache.put(phaseset, bucket);
                return bucket;
              }))
          .orElse(Stream.builder())
          .add(call);
    });
    return partition.build()
        .map(bucket -> bucket.build().collect(Collectors.toList()))
        .collect(Collectors.toList());
  }
}

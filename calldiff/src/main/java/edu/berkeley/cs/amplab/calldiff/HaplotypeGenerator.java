/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

    private static boolean equals(char lhs, char rhs) {
      switch (lhs) {
        case 'A':
          switch (rhs) {
            case 'A':
            case 'R':
            case 'M':
            case 'W':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
              return true;
            case 'C':
            case 'G':
            case 'T':
            case 'U':
            case 'Y':
            case 'K':
            case 'S':
            case 'B':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'C':
          switch (rhs) {
            case 'C':
            case 'Y':
            case 'M':
            case 'S':
            case 'B':
            case 'H':
            case 'V':
            case 'N':
              return true;
            case 'A':
            case 'G':
            case 'T':
            case 'U':
            case 'R':
            case 'K':
            case 'W':
            case 'D':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'G':
          switch (rhs) {
            case 'G':
            case 'R':
            case 'K':
            case 'S':
            case 'B':
            case 'D':
            case 'V':
            case 'N':
              return true;
            case 'A':
            case 'C':
            case 'T':
            case 'U':
            case 'Y':
            case 'M':
            case 'W':
            case 'H':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'T':
          switch (rhs) {
            case 'T':
            case 'Y':
            case 'K':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'N':
              return true;
            case 'A':
            case 'C':
            case 'G':
            case 'U':
            case 'R':
            case 'M':
            case 'S':
            case 'V':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'U':
          switch (rhs) {
            case 'U':
            case 'Y':
            case 'K':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'N':
              return true;
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'R':
            case 'M':
            case 'S':
            case 'V':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'R':
          switch (rhs) {
            case 'R':
            case 'A':
            case 'G':
              return true;
            case 'C':
            case 'T':
            case 'U':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'Y':
          switch (rhs) {
            case 'Y':
            case 'C':
            case 'T':
            case 'U':
              return true;
            case 'A':
            case 'G':
            case 'R':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'K':
          switch (rhs) {
            case 'K':
            case 'G':
            case 'T':
            case 'U':
              return true;
            case 'A':
            case 'C':
            case 'R':
            case 'Y':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'M':
          switch (rhs) {
            case 'M':
            case 'A':
            case 'C':
              return true;
            case 'G':
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'S':
          switch (rhs) {
            case 'S':
            case 'C':
            case 'G':
              return true;
            case 'A':
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'W':
          switch (rhs) {
            case 'W':
            case 'A':
            case 'T':
            case 'U':
              return true;
            case 'C':
            case 'G':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'B':
          switch (rhs) {
            case 'B':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
              return true;
            case 'A':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'D':
          switch (rhs) {
            case 'D':
            case 'A':
            case 'G':
            case 'T':
            case 'U':
              return true;
            case 'C':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'H':
          switch (rhs) {
            case 'H':
            case 'A':
            case 'C':
            case 'T':
            case 'U':
              return true;
            case 'G':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'V':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'V':
          switch (rhs) {
            case 'V':
            case 'A':
            case 'C':
            case 'G':
              return true;
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'N':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'N':
          switch (rhs) {
            case 'N':
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
              return true;
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'X':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case 'X':
          switch (rhs) {
            case 'X':
              return true;
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case '-':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        case '-':
          switch (rhs) {
            case '-':
              return true;
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
              return false;
            default:
              throw new IllegalStateException(Character.toString(rhs));
          }
        default:
          throw new IllegalStateException(Character.toString(lhs));
      }
    }

    private static boolean equals(String lhs, String rhs) {
      int length = lhs.length();
      if (length == rhs.length()) {
        for (int i = 0; i < length; ++i) {
          if (!equals(lhs.charAt(i), rhs.charAt(i))) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    private static <X> X getProperty(List<Call> calls, Function<? super Call, ? extends X> accessor,
        String message) {
      List<X> list = calls.stream().map(accessor).distinct().collect(Collectors.toList());
      Preconditions.checkState(1 == list.size(), message);
      return Iterables.getOnlyElement(list);
    }
    private final NavigableMap<Integer, Integer> offsets;

    private final String string;

    private Haplotype(String string, NavigableMap<Integer, Integer> offsets) {
      this.string = string;
      this.offsets = offsets;
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
                  .map(CaseNormalizer::normalizeCase)
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
    public boolean equals(Object obj) {
      return null != obj && Haplotype.class != obj.getClass()
          && equals(toString(), ((Haplotype) obj).toString());
    }

    @Override
    public int hashCode() {
      return toString().hashCode();
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
    return haplotypes.collect(Collectors.toSet())
        .stream()
        .map(Object::toString)
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

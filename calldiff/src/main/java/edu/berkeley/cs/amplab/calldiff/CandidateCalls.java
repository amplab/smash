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

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.RandomAccess;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A {@code CandidateCalls} object is a 5-tuple containing the following information:
 *
 * <ul>
 *   <li>The contig name that the calls are taken from</li>
 *   <li>The start position of the window (the infimum of all the call positions on both the left
 *       and right hand side)</li>
 *   <li>The end position of the window (the supremum of all the call ends on both the left and
 *       right hand side)</li>
 *   <li>The calls on the window from the left hand side</li>
 *   <li>The calls on the window from the right hand side</li>
 * </ul>
 */
public class CandidateCalls {

  private static final HashCodeAndEquals<CandidateCalls> HASH_CODE_AND_EQUALS =
      HashCodeAndEquals.create(CandidateCalls.class,
          CandidateCalls::contig,
          CandidateCalls::start,
          CandidateCalls::end,
          CandidateCalls::lhs,
          CandidateCalls::rhs);

  public static CandidateCalls create(
      String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
    return new CandidateCalls(contig, start, end, lhs, rhs);
  }

  public static Stream<CandidateCalls> createCandidates(Window window) {
    return window.isTooLarge()
        ? Stream.empty()
        : BimonotonicAStarSearcher.<List<Call>, List<Call>, CandidateCalls>builder()
            .setBiFunction(
                (lhs, rhs) -> create(window.contig(), window.start(), window.end(), lhs, rhs))
            .setComparator(Comparator.comparing(CandidateCalls::size).reversed())
            .build()
            .search(nonOverlappingSubsets(window.lhs()), nonOverlappingSubsets(window.rhs()));
  }

  private static <L extends List<? extends Call> & RandomAccess> ArrayList<List<Call>>
      nonOverlappingSubsets(L calls) {
    ArrayList<List<Call>> list = Sets.powerSet(Sets.newLinkedHashSet(calls))
        .stream()
        .filter(
            set -> {
              for (Call lhs : set) {
                for (Call rhs : set) {
                  if (lhs != rhs && lhs.overlaps(rhs)) {
                    return false;
                  }
                }
              }
              return true;
            })
        .map(Lists::newArrayList)
        .collect(Collectors.toCollection(ArrayList::new));
    Collections.sort(
        list,
        Comparator.comparing((Function<Collection<?>, Integer>) Collection::size).reversed());
    return list;
  }

  private final String contig;
  private final List<Call> lhs, rhs;
  private final int start, end;

  private CandidateCalls(String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
    this.contig = contig;
    this.start = start;
    this.end = end;
    this.lhs = lhs;
    this.rhs = rhs;
  }

  public String contig() {
    return contig;
  }

  public int end() {
    return end;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  public boolean generatesSameSetOfHaplotypes(FastaReader.FastaFile reference) {
    String contig = contig();
    int start = start(), end = end();
    return Objects.equals(
        HaplotypeGenerator.generateHaplotypes(reference, contig, lhs(), start, end),
        HaplotypeGenerator.generateHaplotypes(reference, contig, rhs(), start, end));
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public List<Call> lhs() {
    return lhs;
  }

  public List<Call> rhs() {
    return rhs;
  }

  int size() {
    return lhs().size() + rhs().size();
  }

  public int start() {
    return start;
  }
}

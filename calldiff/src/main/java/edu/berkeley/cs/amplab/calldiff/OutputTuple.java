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

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Contains the output results for a single genomic window. It has four fields:
 *
 * <ul>
 *   <li>The calls from the left hand side that had equivalents on the right hand side</li>
 *   <li>The calls from the left hand side that did not have equivalents on the right hand side</li>
 *   <li>The calls from the right hand side that had equivalents on the left hand side</li>
 *   <li>The calls from the right hand side that did not have equivalents on the left hand side</li>
 * </ul>
 */
public class OutputTuple {

  public static class Builder {

    private static <X> List<X> toList(Stream.Builder<? extends X> stream) {
      return stream.build().collect(Collectors.toList());
    }

    private final Stream.Builder<Call>
        matchingLhs = Stream.builder(),
        matchingRhs = Stream.builder(),
        notMatchingLhs = Stream.builder(),
        notMatchingRhs = Stream.builder();
    private final Window window;

    private Builder(Window window) {
      this.window = window;
    }

    private <X> Builder addAll(Iterable<? extends X> iterable, Stream.Builder<X> stream) {
      for (X object : iterable) {
        stream.add(object);
      }
      return this;
    }

    public Builder addMatchingLhs(Iterable<Call> calls) {
      return addAll(calls, matchingLhs);
    }

    public Builder addMatchingRhs(Iterable<Call> calls) {
      return addAll(calls, matchingRhs);
    }

    public Builder addNotMatchingLhs(Iterable<Call> calls) {
      return addAll(calls, notMatchingLhs);
    }

    public Builder addNotMatchingRhs(Iterable<Call> calls) {
      return addAll(calls, notMatchingRhs);
    }

    public OutputTuple build() {
      return new OutputTuple(
          window,
          toList(matchingLhs),
          toList(matchingRhs),
          toList(notMatchingLhs),
          toList(notMatchingRhs));
    }
  }

  private static final
      HashCodeAndEquals<OutputTuple> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
          OutputTuple.class,
          OutputTuple::window,
          OutputTuple::matchingLhs,
          OutputTuple::matchingRhs,
          OutputTuple::notMatchingLhs,
          OutputTuple::notMatchingRhs);

  public static Builder builder(Window window) {
    return new Builder(window);
  }

  public static Stream<OutputTuple>
      calldiff(FastaReader.FastaFile reference, Stream<Call> lhs, Stream<Call> rhs) {
    return Window.partition(lhs, rhs)
        .map(window -> window.createOutputTuple(window.candidates()
            .filter(candidates -> candidates.generatesSameSetOfHaplotypes(reference))
            .findFirst()));
  }

  private final List<Call> matchingLhs, matchingRhs, notMatchingLhs, notMatchingRhs;
  private final Window window;

  private OutputTuple(
      Window window,
      List<Call> matchingLhs,
      List<Call> matchingRhs,
      List<Call> notMatchingLhs,
      List<Call> notMatchingRhs) {
    this.window = window;
    this.matchingLhs = matchingLhs;
    this.matchingRhs = matchingRhs;
    this.notMatchingLhs = notMatchingLhs;
    this.notMatchingRhs = notMatchingRhs;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public List<Call> matchingLhs() {
    return matchingLhs;
  }

  public List<Call> matchingRhs() {
    return matchingRhs;
  }

  public List<Call> notMatchingLhs() {
    return notMatchingLhs;
  }

  public List<Call> notMatchingRhs() {
    return notMatchingRhs;
  }

  public Window window() {
    return window;
  }
}

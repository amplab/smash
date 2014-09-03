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

import com.google.common.collect.ImmutableMap;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Stream;

/**
 * A {@link Collector} that knows how to create a {@link Map} that remembers the index of each
 * object in the {@link Stream} that it encounters.
 */
public class Indexer<X> implements Collector<X, Indexer<X>, Map<X, Integer>> {

  public static <X> Indexer<X> create() {
    return new Indexer<>();
  }

  private final ImmutableMap.Builder<X, Integer> builder = ImmutableMap.builder();
  private int index = 0;

  private Indexer() {}

  @Override public BiConsumer<Indexer<X>, X> accumulator() {
    return (indexer, object) -> indexer.builder.put(object, index++);
  }

  @Override public Set<Collector.Characteristics> characteristics() {
    return Collections.emptySet();
  }

  @Override public BinaryOperator<Indexer<X>> combiner() {
    return (lhs, rhs) -> {
      throw new UnsupportedOperationException();
    };
  }

  @Override public Function<Indexer<X>, Map<X, Integer>> finisher() {
    return indexer -> indexer.builder.build();
  }

  @Override public Supplier<Indexer<X>> supplier() {
    return Indexer::create;
  }
}

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

import static java.util.Spliterator.DISTINCT;
import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.NONNULL;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Queue;
import java.util.Set;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * An implementation of the
 * <a href="http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm">Bron-Kerbosch</a>
 * algorithm used for finding maximal independent sets in the
 * <a href="http://en.wikipedia.org/wiki/Interval_graph">interval graph</a> described by a set of
 * variant calls.
 */
public class BronKerbosch<V> {

  private static class DepthFirstSearch<X> {

    static class Builder<X> {

      private Function<? super X, ? extends Stream<? extends X>> neighbors;

      DepthFirstSearch<X> build() {
        return new DepthFirstSearch<>(neighbors);
      }

      Builder<X> setNeighborsFunction(
          Function<? super X, ? extends Stream<? extends X>> neighbors) {
        this.neighbors = neighbors;
        return this;
      }
    }

    static <X> Builder<X> builder() {
      return new Builder<>();
    }

    private final Function<? super X, ? extends Stream<? extends X>> neighbors;

    private DepthFirstSearch(Function<? super X, ? extends Stream<? extends X>> neighbors) {
      this.neighbors = neighbors;
    }

    Stream<X> search(X initial) {
      return stream(
          new AbstractIterator<X>() {

            private final Queue<X> stack = Collections.asLifoQueue(new ArrayDeque<>());

            {
              stack.offer(initial);
            }

            @Override protected X computeNext() {
              X next = stack.poll();
              if (null != next) {
                neighbors.apply(next).forEach(stack::offer);
                return next;
              }
              return endOfData();
            }
          });
    }
  }

  public static <V> Stream<Set<V>> search(SetMultimap<V, V> graph) {
    class Frame {

      private Set<V> r, p, x;

      Frame(Set<V> r, Set<V> p, Set<V> x) {
        this.r = r;
        this.p = p;
        this.x = x;
      }

      boolean isMaximal() {
        return p.isEmpty() && x.isEmpty();
      }

      Stream<Frame> neighbors() {
        return isMaximal()
            ? Stream.empty()
            : stream(
                new AbstractIterator<Frame>() {

                  private final Iterator<V> iterator = p.iterator();

                  @Override protected Frame computeNext() {
                    if (iterator.hasNext()) {
                      V vertex = iterator.next();
                      Set<V> neighbors = graph.get(vertex);
                      Frame next = new Frame(
                          copy(Stream.concat(r.stream(), Stream.of(vertex))),
                          copy(Sets.intersection(p, neighbors).stream()),
                          copy(Sets.intersection(x, neighbors).stream()));
                      iterator.remove();
                      x.add(vertex);
                      return next;
                    }
                    return endOfData();
                  }

                  private <X> Set<X> copy(Stream<? extends X> stream) {
                    return stream.collect(Collectors.toCollection(LinkedHashSet::new));
                  }
                });
      }

      Set<V> r() {
        return r;
      }
    }
    return DepthFirstSearch.<Frame>builder()
        .setNeighborsFunction(Frame::neighbors)
        .build()
        .search(new Frame(
            new LinkedHashSet<>(),
            Sets.newLinkedHashSet(graph.keySet()),
            new LinkedHashSet<>()))
        .filter(Frame::isMaximal)
        .map(Frame::r);
  }

  private static <X> Stream<X> stream(Iterator<X> iterator) {
    return StreamSupport.stream(
        Spliterators.spliteratorUnknownSize(iterator, DISTINCT | IMMUTABLE | NONNULL),
        false);
  }
}

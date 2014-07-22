package edu.berkeley.cs.amplab.vardiff;

import static java.util.Spliterator.DISTINCT;
import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.NONNULL;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Optional;
import java.util.Queue;
import java.util.Set;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class ModifiedBronKerbosch<V> {

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

  private static <X> Stream<X> stream(Iterator<X> iterator) {
    return StreamSupport.stream(
        Spliterators.spliteratorUnknownSize(iterator, DISTINCT | IMMUTABLE | NONNULL),
        false);
  }

  public static <V> Stream<Set<V>> search(SetMultimap<V, V> graph) {
    class Frame {

      private Set<V> r, p, x;

      Frame(Set<V> r, Set<V> p, Set<V> x) {
        this.r = r;
        this.p = p;
        this.x = x;
      }

      Stream<Frame> neighbors() {
        return p.isEmpty() && x.isEmpty()
            ? Stream.empty()
            : stream(
                new AbstractIterator<Frame>() {
                  @Override protected Frame computeNext() {
                    Optional<Frame> next = Optional
                        .ofNullable(Iterators.getNext(p.iterator(), null))
                        .map(
                            vertex -> {
                              Set<V>
                                  neighbors = graph.get(vertex),
                                  singleton = Collections.singleton(vertex);
                              Frame frame = new Frame(
                                  Sets.union(r, singleton),
                                  Sets.intersection(p, neighbors),
                                  Sets.intersection(x, neighbors));
                              p = Sets.difference(p, singleton);
                              x = Sets.union(x, singleton);
                              return frame;
                            });
                    return next.isPresent() ? next.get() : endOfData();
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
        .map(Frame::r);
  }
}

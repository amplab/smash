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

package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DiffStats {

  public static class Builder implements Collector<OutputTuple, Builder, DiffStats> {

    private static final Set<Collector.Characteristics>
        CHARACTERISTICS = Collections.singleton(Collector.Characteristics.UNORDERED);

    private static void
        addAll(ImmutableMultiset.Builder<Call.Type> multiset, Collection<? extends Call> calls) {
      multiset.addAll(calls.stream().map(Call::type).iterator());
    }

    private final ImmutableMultiset.Builder<Call.Type>
        matchingLhs = ImmutableMultiset.builder(),
        matchingRhs = ImmutableMultiset.builder(),
        notMatchingLhs = ImmutableMultiset.builder(),
        notMatchingRhs = ImmutableMultiset.builder();
    private final ImmutableList.Builder<Window>
        unprocessedWindows = ImmutableList.builder();

    private Builder() {}

    @Override public BiConsumer<Builder, OutputTuple> accumulator() {
      return (builder, tuple) -> {
        addAll(builder.matchingLhs, tuple.matchingLhs());
        addAll(builder.matchingRhs, tuple.matchingRhs());
        addAll(builder.notMatchingLhs, tuple.notMatchingLhs());
        addAll(builder.notMatchingRhs, tuple.notMatchingRhs());
        Window window = tuple.window();
        if (window.isTooLarge()) {
          unprocessedWindows.add(window);
        }
      };
    }

    @Override public Set<Collector.Characteristics> characteristics() {
      return CHARACTERISTICS;
    }

    @Override public BinaryOperator<Builder> combiner() {
      return (lhs, rhs) -> {
        throw new UnsupportedOperationException();
      };
    }

    @Override public Function<Builder, DiffStats> finisher() {
      return builder -> new DiffStats(
          builder.matchingLhs.build(),
          builder.matchingRhs.build(),
          builder.notMatchingLhs.build(),
          builder.notMatchingRhs.build(),
          builder.unprocessedWindows.build());
    }

    @Override public Supplier<Builder> supplier() {
      return DiffStats::builder;
    }
  }

  private static final HashCodeAndEquals<DiffStats> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
      DiffStats.class,
      DiffStats::matchingLhs,
      DiffStats::matchingRhs,
      DiffStats::notMatchingLhs,
      DiffStats::notMatchingRhs);

  private static final String NEWLINE = String.format("%n");

  public static Builder builder() {
    return new Builder();
  }
  private final Multiset<Call.Type> matchingLhs, matchingRhs, notMatchingLhs, notMatchingRhs;

  private final List<Window> unprocessedWindows;

  private DiffStats(
      Multiset<Call.Type> matchingLhs,
      Multiset<Call.Type> matchingRhs,
      Multiset<Call.Type> notMatchingLhs,
      Multiset<Call.Type> notMatchingRhs,
      List<Window> unprocessedWindows) {
    this.matchingLhs = matchingLhs;
    this.matchingRhs = matchingRhs;
    this.notMatchingLhs = notMatchingLhs;
    this.notMatchingRhs = notMatchingRhs;
    this.unprocessedWindows = unprocessedWindows;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public Multiset<Call.Type> matchingLhs() {
    return matchingLhs;
  }

  public Multiset<Call.Type> matchingRhs() {
    return matchingRhs;
  }

  public Multiset<Call.Type> notMatchingLhs() {
    return notMatchingLhs;
  }

  public Multiset<Call.Type> notMatchingRhs() {
    return notMatchingRhs;
  }

  @Override
  public String toString() {
    return Stream
        .of(
            Stream
                .concat(
                    Stream.of(Stream.of(
                        "TYPE",
                        "MATCHING LHS",
                        "MATCHING RHS",
                        "NOT MATCHING LHS",
                        "NOT MATCHING RHS")),
                    Stream.of(Call.Type.values()).map(type -> Stream.of(
                        type,
                        matchingLhs.count(type),
                        matchingRhs.count(type),
                        notMatchingLhs.count(type),
                        notMatchingRhs.count(type))))
                .map(stream -> stream.map(Object::toString).collect(Collectors.joining("\t")))
                .collect(Collectors.joining(NEWLINE)),
            String.format(
                "Unprocessed windows: %s",
                unprocessedWindows.stream()
                    .map(Object::toString)
                    .collect(Collectors.joining(", "))))
        .collect(Collectors.joining(NEWLINE));
  }

  public List<Window> unprocessedWindows() {
    return unprocessedWindows;
  }
}

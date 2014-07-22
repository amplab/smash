package edu.berkeley.cs.amplab.calldiff;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

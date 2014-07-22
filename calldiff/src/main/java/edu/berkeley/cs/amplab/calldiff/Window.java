package edu.berkeley.cs.amplab.calldiff;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Sets;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Window {

  public static class Builder {

    private final String contig;
    private ArrayList<Call>
        lhs = new ArrayList<>(),
        rhs = new ArrayList<>();
    private int
        start = Integer.MAX_VALUE,
        end = Integer.MIN_VALUE;

    Builder(String contig) {
      this.contig = contig;
    }

    public Builder addLhs(Call call) {
      lhs.add(call);
      return updateStartAndEnd(call);
    }

    public Builder addRhs(Call call) {
      rhs.add(call);
      return updateStartAndEnd(call);
    }

    public Window build() {
      return new Window(contig, start, end, lhs, rhs);
    }

    int end() {
      return end;
    }

    private Builder updateStartAndEnd(Call call) {
      start = Math.min(start, call.position());
      end = Math.max(end, call.end());
      return this;
    }
  }

  private static class CallWithSource {

    static final Comparator<CallWithSource> COMPARATOR = Comparator
        .<CallWithSource, String>comparing(callWithSource -> callWithSource.call().contig())
        .thenComparing(callWithSource -> callWithSource.call().position());

    private final Call call;
    private final Source source;

    CallWithSource(Call call, Source source) {
      this.call = call;
      this.source = source;
    }

    Call call() {
      return call;
    }

    Source source() {
      return source;
    }
  }

  private enum Source {

    LHS,
    RHS;

    Iterator<CallWithSource> iterator(Stream<Call> calls) {
      return calls.map(call -> new CallWithSource(call, this)).iterator();
    }
  }

  private static final HashCodeAndEquals<Window> HASH_CODE_AND_EQUALS =
      HashCodeAndEquals.create(
          Window.class,
          Window::contig,
          Window::start,
          Window::end,
          Window::lhs,
          Window::rhs);

  private static final int MAX_WINDOW_SIZE = 10;

  public static Builder builder(String contig) {
    return new Builder(contig);
  }

  private static <X> Comparator<X> comparator(Collection<? extends X> collection) {
    return Comparator.comparing(collection.stream().collect(Indexer.create())::get);
  }

  public static Window
      create(String contig, int start, int end, ArrayList<Call> lhs, ArrayList<Call> rhs) {
    return new Window(contig, start, end, lhs, rhs);
  }

  public static Stream<Window> partition(final Stream<Call> lhs, final Stream<Call> rhs) {
    return StreamSupport.stream(
        Spliterators.spliteratorUnknownSize(
            new AbstractIterator<Window>() {

              private final PeekingIterator<CallWithSource> iterator = Iterators.peekingIterator(
                  Iterators.mergeSorted(
                      Arrays.asList(Source.LHS.iterator(lhs), Source.RHS.iterator(rhs)),
                      CallWithSource.COMPARATOR));

              private void addToWindow(Window.Builder window, CallWithSource next, Call call) {
                switch (next.source()) {
                  case LHS:
                    window.addLhs(call);
                    break;
                  case RHS:
                    window.addRhs(call);
                    break;
                  default:
                    throw new IllegalStateException();
                }
              }

              @Override protected Window computeNext() {
                if (iterator.hasNext()) {
                  CallWithSource next = iterator.next();
                  Call call = next.call();
                  String contig = call.contig();
                  Window.Builder window = Window.builder(contig);
                  for (addToWindow(window, next, call);
                      iterator.hasNext()
                          && Objects.equals(
                              contig,
                              (call = (next = iterator.peek()).call()).contig())
                          && call.position() < window.end();) {
                    addToWindow(window, iterator.next(), call);
                  }
                  return window.build();
                }
                return endOfData();
              }
            },
            Spliterator.DISTINCT | Spliterator.IMMUTABLE | Spliterator.NONNULL),
        false);
  }

  private static <X> Set<X> set(Iterable<? extends X> iterable) {
    Set<X> set = new HashSet<>();
    for (X object : iterable) {
      set.add(object);
    }
    return set;
  }

  private static <X> List<X> sort(Collection<? extends X> collection,
      Comparator<? super X> comparator) {
    List<X> list = new ArrayList<>(collection.size());
    list.addAll(collection);
    Collections.sort(list, comparator);
    return list;
  }

  private final String contig;
  private final ArrayList<Call> lhs, rhs;
  private final int start, end;

  private Window(String contig, int start, int end, ArrayList<Call> lhs, ArrayList<Call> rhs) {
    this.contig = contig;
    this.start = start;
    this.end = end;
    this.lhs = lhs;
    this.rhs = rhs;
  }

  public Stream<CandidateCalls> candidates() {
    return CandidateCalls.createCandidates(this);
  }

  public String contig() {
    return contig;
  }

  public OutputTuple createOutputTuple(Optional<CandidateCalls> optional) {
    OutputTuple.Builder tuple = OutputTuple.builder(this);
    if (optional.isPresent()) {
      CandidateCalls calls = optional.get();
      List<Call>
          windowLhs = lhs(),
          windowRhs = rhs();
      Comparator<Call>
          lhsComparator = comparator(windowLhs),
          rhsComparator = comparator(windowRhs);
      Set<Call>
          windowLhsSet = set(windowLhs),
          windowRhsSet = set(windowRhs),
          callsLhsSet = set(calls.lhs()),
          callsRhsSet = set(calls.rhs());
      tuple
          .addMatchingLhs(sort(Sets.intersection(windowLhsSet, callsLhsSet), lhsComparator))
          .addMatchingRhs(sort(Sets.intersection(windowRhsSet, callsRhsSet), rhsComparator))
          .addNotMatchingLhs(sort(Sets.difference(windowLhsSet, callsLhsSet), lhsComparator))
          .addNotMatchingRhs(sort(Sets.difference(windowRhsSet, callsRhsSet), rhsComparator));
    }
    return tuple.build();
  }

  public int end() {
    return end;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public boolean isTooLarge() {
    return MAX_WINDOW_SIZE < size();
  }

  public ArrayList<Call> lhs() {
    return lhs;
  }

  public ArrayList<Call> rhs() {
    return rhs;
  }

  public int size() {
    return Math.max(lhs().size(), rhs.size());
  }

  public int start() {
    return start;
  }

  @Override
  public String toString() {
    return Stream
        .of(
            Stream.of("contig", contig()),
            Stream.of("start", start()),
            Stream.of("end", end()),
            Stream.of("lhs", lhs()),
            Stream.of("rhs", rhs()))
        .map(stream -> stream.map(Object::toString).collect(Collectors.joining("=")))
        .collect(Collectors.joining(", "));
  }
}

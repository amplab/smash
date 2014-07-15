package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Window {

  public static class Builder {

    private final String contig;
    private ImmutableList.Builder<Call>
        lhs = ImmutableList.builder(),
        rhs = ImmutableList.builder();
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
      return new Window(contig, start, end, lhs.build(), rhs.build());
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

  public static Builder builder(String contig) {
    return new Builder(contig);
  }

  public static Window create(String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
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

              private Window.Builder addToWindow(Window.Builder window, CallWithSource next, Call call) {
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
                return window;
              }

              @Override protected Window computeNext() {
                if (iterator.hasNext()) {
                  CallWithSource next = iterator.next();
                  Call call = next.call();
                  String contig = call.contig();
                  Window.Builder window = addToWindow(Window.builder(contig), next, call);
                  while (iterator.hasNext()) {
                    next = iterator.peek();
                    call = next.call();
                    if (!Objects.equals(contig, call.contig()) || window.end() <= call.position()) {
                      break;
                    }
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

  private final String contig;
  private final List<Call> lhs, rhs;
  private final int start, end;

  private Window(String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
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

  public List<Call> lhs() {
    return lhs;
  }

  public List<Call> rhs() {
    return rhs;
  }

  public int start() {
    return start;
  }
}

package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class GenomePartitioner {

  private static class CallWithSource implements Comparable<CallWithSource> {

    private static final Comparator<CallWithSource>
        COMPARATOR = Comparator.comparing(CallWithSource::call);

    private final Call call;
    private final Source source;

    CallWithSource(Call call, Source source) {
      this.call = call;
      this.source = source;
    }

    Call call() {
      return call;
    }

    @Override public int compareTo(CallWithSource rhs) {
      return COMPARATOR.compare(this, rhs);
    }

    Source source() {
      return source;
    }
  }

  private enum Source {

    LHS,
    RHS;

    CallWithSource create(Call call) {
      return new CallWithSource(call, this);
    }
  }

  public static class Window {

    static class Builder {

      private int beginning = Integer.MAX_VALUE;
      private final String contig;
      private int end = Integer.MIN_VALUE;
      private final ImmutableList.Builder<Call> lhs = ImmutableList.builder();
      private final ImmutableList.Builder<Call> rhs = ImmutableList.builder();

      private Builder(String contig) {
        this.contig = contig;
      }

      Builder add(CallWithSource callWithSource) {
        Call call = callWithSource.call();
        switch (callWithSource.source()) {
          case LHS:
            lhs.add(call);
            break;
          case RHS:
            rhs.add(call);
            break;
          default:
            throw new IllegalStateException();
        }
        int position = call.position();
        beginning = Math.min(beginning, position);
        end = Math.max(end, position + call.reference().length());
        return this;
      }

      Window build() {
        return new Window(contig, lhs.build(), rhs.build(), beginning, end);
      }

      Builder setBeginning(int beginning) {
        this.beginning = beginning;
        return this;
      }

      Builder setEnd(int end) {
        this.end = end;
        return this;
      }

      Builder setLhs(List<Call> lhs) {
        this.lhs.addAll(lhs);
        return this;
      }

      Builder setRhs(List<Call> rhs) {
        this.rhs.addAll(rhs);
        return this;
      }
    }

    private static final HashCodeAndEquals<Window> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
        Window.class,
        Window::contig,
        Window::lhs,
        Window::rhs,
        Window::beginning,
        Window::end);

    static Builder builder(String contig) {
      return new Builder(contig);
    }

    private final int beginning;
    private final String contig;
    private final int end;
    private final List<Call> lhs;
    private final List<Call> rhs;

    private Window(String contig, List<Call> lhs, List<Call> rhs, int beginning, int end) {
      this.contig = contig;
      this.lhs = lhs;
      this.rhs = rhs;
      this.beginning = beginning;
      this.end = end;
    }

    public int beginning() {
      return beginning;
    }

    public String contig() {
      return contig;
    }

    public int end() {
      return end;
    }

    @Override public boolean equals(Object obj) {
      return HASH_CODE_AND_EQUALS.equals(this, obj);
    }

    @Override public int hashCode() {
      return HASH_CODE_AND_EQUALS.hashCode(this);
    }

    public List<Call> lhs() {
      return lhs;
    }

    public List<Call> rhs() {
      return rhs;
    }
  }

  public static Stream<Window> partition(final Stream<Call> lhs, final Stream<Call> rhs) {
    return StreamSupport.stream(
        Spliterators.spliteratorUnknownSize(
            new AbstractIterator<Window>() {

              private final PeekingIterator<CallWithSource> iterator = Iterators.peekingIterator(
                  Iterators.mergeSorted(
                      Stream.of(lhs.map(Source.LHS::create), rhs.map(Source.RHS::create))
                          .map(Stream::iterator)
                          .collect(Collectors.toList()),
                      Comparator.naturalOrder()));

              @Override protected Window computeNext() {
                if (iterator.hasNext()) {
                  CallWithSource next = iterator.next();
                  Call call = next.call();
                  int end = call.position() + call.reference().length();
                  String contig = call.contig();
                  Window.Builder window = Window.builder(contig).add(next);
                  while (iterator.hasNext()) {
                    Call peek = iterator.peek().call();
                    if (Objects.equals(contig, peek.contig()) && peek.position() < end) {
                      end = Math.max(end, (call = (next = iterator.next()).call()).position()
                          + call.reference().length());
                      window.add(next);
                    } else {
                      break;
                    }
                  }
                  return window.build();
                }
                return endOfData();
              }
            },
            Spliterator.DISTINCT | Spliterator.IMMUTABLE | Spliterator.NONNULL),
        false);
  }
}

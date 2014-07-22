package edu.berkeley.cs.amplab.vardiff;

import static com.google.common.collect.Iterators.peekingIterator;
import static java.lang.Integer.MAX_VALUE;
import static java.lang.Integer.MIN_VALUE;
import static java.lang.Integer.max;
import static java.lang.Integer.min;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.Stream;

public class CandidateCallsTest {

  private static final Comparator<Call> COMPARATOR = Comparator.comparing(Call::position);

  private static <X, C extends Collection<? extends X>> C
      assertSorted(C collection, Comparator<? super X> comparator) {
    for (
        PeekingIterator<? extends X> iterator = Iterators.peekingIterator(collection.iterator());
        iterator.hasNext();) {
      X next = iterator.next();
      if (iterator.hasNext()) {
        assertTrue(comparator.compare(next, iterator.peek()) <= 0);
      }
    }
    return collection;
  }

  private static List<Call> assertSorted(List<Call> calls) {
    return assertSorted(calls, COMPARATOR);
  }

  @Test
  public void testCreateCandidates() throws IOException {
    TestReference.reader().read(reference -> {
      String contig = reference.contigs().iterator().next();
      int contigLen = reference.contigLength(contig);
      Random random = new Random();
      for (int i = 1; i < 40; ++i) {
        final int maxCallLen = i;
        for (int j = 1; j < 10; ++j) {
          final int lhsNumCalls = j;
          for (int k = 1; k < 10; ++k) {
            final int rhsNumCalls = k;
              ArrayList<Call>
                  lhs = TestCall.randomCalls(random, contig, contigLen, maxCallLen, lhsNumCalls),
                  rhs = TestCall.randomCalls(random, contig, contigLen, maxCallLen, rhsNumCalls);
              for (
                  PeekingIterator<CandidateCalls> iterator = peekingIterator(CandidateCalls
                      .createCandidates(
                          Stream.concat(lhs.stream(), rhs.stream()).collect(
                              () -> Window.create(contig, MAX_VALUE, MIN_VALUE, lhs, rhs),
                              (window, call) -> Window.create(
                                  window.contig(),
                                  min(window.start(), call.position()),
                                  max(window.end(), call.end()),
                                  window.lhs(),
                                  window.rhs()),
                              (l, r) -> { throw new UnsupportedOperationException(); }))
                      .iterator());
                  iterator.hasNext();) {
                CandidateCalls next = iterator.next();
                for (List<Call> calls :
                    Arrays.asList(assertSorted(next.lhs()), assertSorted(next.rhs()))) {
                  for (Call call1 : calls) {
                    for (Call call2 : calls) {
                      if (call1 != call2 && call1.overlaps(call2)) {
                        fail();
                      }
                    }
                  }
                }
                if (iterator.hasNext()) {
                  assertTrue(iterator.peek().size() <= next.size());
                }
              }
          }
        }
      }
      return null;
    });
  }
}

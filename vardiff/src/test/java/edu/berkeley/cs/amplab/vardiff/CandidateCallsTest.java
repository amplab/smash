package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class CandidateCallsTest {

  @Test
  public void testCreateCandidates() throws IOException {
    for (int i = 1; i < 40; ++i) {
      final int maxCallLength = i;
      for (int j = 1; j < 10; ++j) {
        final int lhsNumberOfCalls = j;
        for (int k = 1; k < 10; ++k) {
          final int rhsNumberOfCalls = k;
          TestReference.reader().read(reference -> {
            String contig = reference.contigs().iterator().next();
            int contigLength = reference.contigLength(contig);
            Random random = new Random();
            for (
                PeekingIterator<CandidateCalls> iterator = Iterators.peekingIterator(CandidateCalls
                    .createCandidates(
                        TestCall.randomCalls(
                            random, contig, contigLength, maxCallLength, lhsNumberOfCalls),
                        TestCall.randomCalls(
                            random, contig, contigLength, maxCallLength, rhsNumberOfCalls))
                    .iterator());
                iterator.hasNext();) {
              CandidateCalls next = iterator.next();
              for (List<Call> calls : Arrays.asList(next.lhs(), next.rhs())) {
                TestCall.checkNonOverlapping(calls);
              }
              if (iterator.hasNext()) {
                assertTrue(iterator.peek().size() <= next.size());
              }
            }
            return null;
          });
        }
      }
    }
  }
}

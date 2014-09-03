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

import edu.berkeley.cs.amplab.calldiff.Call;
import edu.berkeley.cs.amplab.calldiff.CandidateCalls;
import edu.berkeley.cs.amplab.calldiff.Window;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.Stream;

/**
 * Unit test for {@link CandidateCalls}
 */
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
      // Create random instances where the maximum length of calls is maxCallLen
      for (int i = 1; i < 40; ++i) {
        final int maxCallLen = i;
        // Create random instances where there are lhsNumCalls number of calls on the left hand
        // side
        for (int j = 1; j < 10; ++j) {
          final int lhsNumCalls = j;
          // Create random instances where there are rhsNumCalls number of calls on the right hand
          // side
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
                // Assert that the candidate calls on the left and right side are both sorted
                // appropriately
                for (List<Call> calls :
                    Arrays.asList(assertSorted(next.lhs()), assertSorted(next.rhs()))) {
                  // Assert that the no pair of calls from either the left hand side or right hand
                  // side overlap.
                  for (Call call1 : calls) {
                    for (Call call2 : calls) {
                      if (call1 != call2 && call1.overlaps(call2)) {
                        fail();
                      }
                    }
                  }
                }
                // Assert that we consider subsets in monotonically nonincreasing order.
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

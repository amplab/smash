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

import static edu.berkeley.cs.amplab.calldiff.TestCall.randomCalls;
import static org.junit.Assert.fail;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.Call;
import edu.berkeley.cs.amplab.calldiff.Window;

import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Unit test for {@link Window}
 */
public class WindowTest {

  @Test
  public void testPartition() throws IOException {
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
              for (
                  PeekingIterator<Window> iterator = Iterators
                      .peekingIterator(Window
                          .partition(
                              randomCalls(random, contig, contigLen, maxCallLen, lhsNumCalls)
                                  .stream(),
                              randomCalls(random, contig, contigLen, maxCallLen, rhsNumCalls)
                                  .stream())
                          .collect(Collectors.toList())
                      .iterator());
                  iterator.hasNext();) {
                List<Call> window = allCalls(iterator.next());
                // Assert that at least one other call on the window overlaps.
                OUTER: for (Call call1 : window) {
                  for (Call call2 : window) {
                    if (call1 == call2 || call1.overlaps(call2)) {
                      continue OUTER;
                    }
                  }
                  fail();
                }
                // Assert no overlapping calls between adjacent windows.
                if (iterator.hasNext()) {
                  List<Call> peek = allCalls(iterator.peek());
                  for (Call call1 : window) {
                    for (Call call2 : peek) {
                      if (call1.overlaps(call2)) {
                        fail();
                      }
                    }
                  }
                }
              }
          }
        }
      }
      return null;
    });
  }

  private static List<Call> allCalls(Window window) {
    return Stream.concat(window.lhs().stream(), window.rhs().stream())
        .collect(Collectors.toList());
  }
}

package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.fail;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import org.junit.Test;

import edu.berkeley.cs.amplab.vardiff.GenomePartitioner.Window;

import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GenomePartitionerTest {

  @Test
  public void testPartition() throws IOException {
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
                PeekingIterator<Window> iterator = Iterators
                    .peekingIterator(
                        GenomePartitioner
                            .partition(
                                TestCall
                                    .randomCalls(
                                        random,
                                        contig,
                                        contigLength,
                                        maxCallLength,
                                        lhsNumberOfCalls)
                                    .stream(),
                                TestCall
                                    .randomCalls(
                                        random,
                                        contig,
                                        contigLength,
                                        maxCallLength,
                                        rhsNumberOfCalls)
                                    .stream())
                            .collect(Collectors.toList())
                    .iterator());
                iterator.hasNext();) {
              List<Call> window = allCalls(iterator.next());
              OUTER: for (Call call1 : window) {
                for (Call call2 : window) {
                  if (call1 == call2 || TestCall.overlaps(call1, call2)) {
                    continue OUTER;
                  }
                }
                fail();
              }
              if (iterator.hasNext()) {
                TestCall.checkNonOverlapping(window, allCalls(iterator.peek()));
              }
            }
            return null;
          });
        }
      }
    }
  }

  private static List<Call> allCalls(Window window) {
    return Stream.concat(window.lhs().stream(), window.rhs().stream())
        .collect(Collectors.toList());
  }
}

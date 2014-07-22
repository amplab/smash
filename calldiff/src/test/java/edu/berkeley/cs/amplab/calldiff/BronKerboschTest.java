package edu.berkeley.cs.amplab.calldiff;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSetMultimap;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.BronKerbosch;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.stream.Collectors;

public class BronKerboschTest {

  @Test
  public void testSearch() {
    Set<Set<Integer>> expected = new HashSet<>();
    expected.add(ImmutableSet.of(2, 3));
    expected.add(ImmutableSet.of(3, 4));
    expected.add(ImmutableSet.of(4, 5));
    expected.add(ImmutableSet.of(4, 6));
    expected.add(ImmutableSet.of(1, 2, 5));
    for (
        Iterator<Set<Integer>> iterator = BronKerbosch
            .search(ImmutableSetMultimap.<Integer, Integer>builder()
                .put(1, 2)
                .put(1, 5)
                .put(2, 1)
                .put(2, 3)
                .put(2, 5)
                .put(3, 2)
                .put(3, 4)
                .put(4, 3)
                .put(4, 5)
                .put(4, 6)
                .put(5, 1)
                .put(5, 2)
                .put(5, 4)
                .put(6, 4)
                .build())
            .collect(Collectors.toSet())
            .iterator();
        iterator.hasNext();) {
      assertTrue(expected.remove(iterator.next()));
    }
    assertTrue(expected.isEmpty());
  }
}

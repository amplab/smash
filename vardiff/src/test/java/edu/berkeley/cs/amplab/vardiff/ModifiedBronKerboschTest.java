package edu.berkeley.cs.amplab.vardiff;

import static edu.berkeley.cs.amplab.vardiff.ModifiedBronKerbosch.search;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.SetMultimap;

import org.junit.Test;

import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.stream.Collectors;

public class ModifiedBronKerboschTest {

  @Test
  public void testSearch() {
    SetMultimap<Integer, Integer> graph = ImmutableSetMultimap.<Integer, Integer>builder()
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
        .build();
    Set<Set<Integer>> expected = new HashSet<>();
    expected.add(Collections.emptySet());
    graph.keySet()
        .stream()
        .map(Collections::singleton)
        .forEach(expected::add);
    graph.entries()
        .stream()
        .filter(entry -> entry.getKey() < entry.getValue())
        .map(entry -> ImmutableSet.of(entry.getKey(), entry.getValue()))
        .forEach(expected::add);
    expected.add(ImmutableSet.of(1, 2, 5));
    for (
        Iterator<Set<Integer>> iterator = search(graph).collect(Collectors.toSet()).iterator();
        iterator.hasNext();) {
      assertTrue(expected.remove(iterator.next()));
    }
    assertTrue(expected.isEmpty());
  }
}

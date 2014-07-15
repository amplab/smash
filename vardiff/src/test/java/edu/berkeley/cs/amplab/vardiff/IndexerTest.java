package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ImmutableMap;

import org.junit.Test;

import java.util.stream.Stream;

public class IndexerTest {

  @Test
  public void testIndexer() {
    assertEquals(
        ImmutableMap.of("foo", 0, "bar", 1, "baz", 2, "fizz", 3, "buzz", 4),
        Stream.of("foo", "bar", "baz", "fizz", "buzz").collect(Indexer.create()));
  }
}

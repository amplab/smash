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

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ImmutableMap;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.Indexer;

import java.util.stream.Stream;

/**
 * Unit test for {@link Indexer}
 */
public class IndexerTest {

  @Test
  public void testIndexer() {
    assertEquals(
        ImmutableMap.of("foo", 0, "bar", 1, "baz", 2, "fizz", 3, "buzz", 4),
        Stream.of("foo", "bar", "baz", "fizz", "buzz").collect(Indexer.create()));
  }
}

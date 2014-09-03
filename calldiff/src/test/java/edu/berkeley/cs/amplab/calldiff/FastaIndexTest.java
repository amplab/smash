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

import com.google.common.collect.ImmutableSortedSet;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.FastaIndex;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 * Unit test for {@link FastaIndex}
 */
public class FastaIndexTest {

  @Test
  public void testFastaIndex() throws IOException {
    FastaIndex expected =
        FastaIndex.createFromEntries(ImmutableSortedSet.<FastaIndex.Entry>naturalOrder()
            .add(FastaIndex.Entry.create("1", 249250621L, 52L, 60, 61))
            .add(FastaIndex.Entry.create("2", 243199373L, 253404903L, 60, 61))
            .add(FastaIndex.Entry.create("3", 198022430L, 500657651L, 60, 61))
            .add(FastaIndex.Entry.create("4", 191154276L, 701980507L, 60, 61))
            .add(FastaIndex.Entry.create("5", 180915260L, 896320740L, 60, 61))
            .add(FastaIndex.Entry.create("6", 171115067L, 1080251307L, 60, 61))
            .add(FastaIndex.Entry.create("7", 159138663L, 1254218344L, 60, 61))
            .add(FastaIndex.Entry.create("8", 146364022L, 1416009371L, 60, 61))
            .add(FastaIndex.Entry.create("9", 141213431L, 1564812846L, 60, 61))
            .add(FastaIndex.Entry.create("10", 135534747L, 1708379889L, 60, 61))
            .add(FastaIndex.Entry.create("11", 135006516L, 1846173603L, 60, 61))
            .add(FastaIndex.Entry.create("12", 133851895L, 1983430282L, 60, 61))
            .add(FastaIndex.Entry.create("13", 115169878L, 2119513096L, 60, 61))
            .add(FastaIndex.Entry.create("14", 107349540L, 2236602526L, 60, 61))
            .add(FastaIndex.Entry.create("15", 102531392L, 2345741279L, 60, 61))
            .add(FastaIndex.Entry.create("16", 90354753L, 2449981581L, 60, 61))
            .add(FastaIndex.Entry.create("17", 81195210L, 2541842300L, 60, 61))
            .add(FastaIndex.Entry.create("18", 78077248L, 2624390817L, 60, 61))
            .add(FastaIndex.Entry.create("19", 59128983L, 2703769406L, 60, 61))
            .add(FastaIndex.Entry.create("20", 63025520L, 2763883926L, 60, 61))
            .add(FastaIndex.Entry.create("21", 48129895L, 2827959925L, 60, 61))
            .add(FastaIndex.Entry.create("22", 51304566L, 2876892038L, 60, 61))
            .add(FastaIndex.Entry.create("X", 155270560L, 2929051733L, 60, 61))
            .add(FastaIndex.Entry.create("Y", 59373566L, 3086910193L, 60, 61))
            .add(FastaIndex.Entry.create("MT", 16569L, 3147273397L, 70, 71))
            .build());
    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
    expected.write(buffer);
    assertEquals(expected, FastaIndex.read(new ByteArrayInputStream(buffer.toByteArray())));
  }
}

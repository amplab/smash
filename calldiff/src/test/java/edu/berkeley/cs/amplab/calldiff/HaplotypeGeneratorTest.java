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

import static edu.berkeley.cs.amplab.calldiff.Call.Phaseset.DEFAULT;
import static edu.berkeley.cs.amplab.calldiff.HaplotypeGenerator.generateHaplotypes;
import static edu.berkeley.cs.amplab.calldiff.HaplotypeGenerator.partitionByPhaseset;
import static edu.berkeley.cs.amplab.calldiff.TestCall.create;
import static edu.berkeley.cs.amplab.calldiff.TestCall.randomCalls;
import static edu.berkeley.cs.amplab.calldiff.TestReference.reference;
import static java.util.Arrays.asList;
import static java.util.Collections.singletonList;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableMap;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.Call;
import edu.berkeley.cs.amplab.calldiff.Call.Phaseset;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Unit test for {@link HaplotypeGenerator}
 */
public class HaplotypeGeneratorTest {

  @Test
  public void testPartitionByPhaseset() throws IOException {
    TestReference.reader().read(reference -> {
      String contig = reference.contigs().iterator().next();
      int contigLength = reference.contigLength(contig);
      Random random = new Random();
      List<Optional<Call.Phaseset>> phasesets = Arrays.asList(
          Optional.empty(),
          Optional.of(DEFAULT),
          Optional.of(Call.Phaseset.create(1)),
          Optional.of(Call.Phaseset.create(2)));
      int numberOfPhasesets = phasesets.size();
      for (List<Call> calls : partitionByPhaseset(
          randomCalls(
              random,
              contig,
              contigLength,
              6,
              20,
              () -> phasesets.get(random.nextInt(numberOfPhasesets))))) {
        Set<Optional<Call.Phaseset>> phasesetsSoFar = new HashSet<>();
        Optional<Optional<Call.Phaseset>> bucket = Optional.empty();
        for (Call call : calls) {
          Optional<Phaseset> phaseset = call.phaseset();
          if (bucket.isPresent()) {
            assertEquals(bucket.get(), phaseset);
          } else {
            bucket = Optional.of(phaseset);
          }
        }
        assertTrue(bucket.isPresent() && phasesetsSoFar.add(bucket.get()));
      }
      return null;
    });
  }

  @Test
  public void testGenerateHaplotypes() {
    String contig = "chr1";
    assertEquals(
        Stream.of("AGGTACTT", "ACGTCCGT", "ACGTACGT", "AGGTCCTT").collect(Collectors.toSet()),
        generateHaplotypes(
            reference(ImmutableMap.of(contig, "ACGTACGT")),
            contig,
            asList(
                create(contig, 1, "C", singletonList("G"), asList(0, 1), DEFAULT),
                create(contig, 4, "A", singletonList("C"), asList(0, 1)),
                create(contig, 6, "G", singletonList("T"), asList(0, 1), DEFAULT)),
            0,
            8));
  }
}

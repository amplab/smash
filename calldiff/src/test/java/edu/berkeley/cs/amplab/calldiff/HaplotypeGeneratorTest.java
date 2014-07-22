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
                create(contig, 2, "C", singletonList("G"), asList(0, 1), DEFAULT),
                create(contig, 5, "A", singletonList("C"), asList(0, 1)),
                create(contig, 7, "G", singletonList("T"), asList(0, 1), DEFAULT)),
            1,
            9));
  }
}

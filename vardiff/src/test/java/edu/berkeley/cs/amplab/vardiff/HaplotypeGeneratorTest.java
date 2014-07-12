package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import edu.berkeley.cs.amplab.vardiff.Call.Phaseset;
import edu.berkeley.cs.amplab.vardiff.HaplotypeGenerator.Haplotype;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.Set;

public class HaplotypeGeneratorTest {

  @Test
  public void testPartitionByPhaseset() throws IOException {
    TestReference.reader().read(reference -> {
      String contig = reference.contigs().iterator().next();
      int contigLength = reference.contigLength(contig);
      Random random = new Random();
      List<Optional<Call.Phaseset>> phasesets = Arrays.asList(
          Optional.empty(),
          Optional.of(Call.Phaseset.DEFAULT),
          Optional.of(Call.Phaseset.create(1)),
          Optional.of(Call.Phaseset.create(2)));
      int numberOfPhasesets = phasesets.size();
      for (List<Call> calls : HaplotypeGenerator.partitionByPhaseset(
          TestCall.randomCalls(
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
  public void testContig_apply() {
    Haplotype.create("ACGTACGT", 5000)
        .apply(Arrays.asList(
            TestCall.create("chr1", 5002, "GT", Collections.singletonList("C"), Arrays.asList(0, 1), Call.Phaseset.DEFAULT),
            TestCall.create("chr1", 5006, "G", Collections.singletonList("A"), Arrays.asList(0, 1), Call.Phaseset.DEFAULT)))
        .map(Object::toString)
        .forEach(System.out::println);
  }
}

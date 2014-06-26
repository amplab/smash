package edu.berkeley.cs.amplab.smash4j;

import static edu.berkeley.cs.amplab.smash4j.TestUtils.getReference;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variant;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variants;
import static org.junit.Assert.assertEquals;

import com.google.common.base.Optional;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.SequenceRescuer.RescuedVariants;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.NavigableMap;

public class SequenceRescuerTest {

  private static File reference;

  private static String multiply(int n, String string) {
    StringBuilder builder = new StringBuilder();
    for (int i = 0; i < n; ++i) {
      builder.append(string);
    }
    return builder.toString();
  }

  @BeforeClass
  public static void setUp() throws IOException {
    reference = getReference();
  }

  private static Optional<SequenceRescuer.RescuedVariants> tryRescue(
      final String contig,
      final int position,
      final NavigableMap<Integer, VariantProto> falseNegatives,
      final NavigableMap<Integer, VariantProto> falsePositives,
      final NavigableMap<Integer, VariantProto> truePositives,
      final int rescueWindowSize) throws Exception {
    return FastaReader.create(reference).read(
        new FastaReader.Callback<Optional<RescuedVariants>>() {
          @Override public Optional<SequenceRescuer.RescuedVariants> read(
              Map<String, Integer> info,
              FastaReader.Callback.FastaFile reference) {
            return SequenceRescuer.builder()
                .setContig(contig)
                .setTruePositives(truePositives)
                .setFalsePositives(falsePositives)
                .setFalseNegatives(falseNegatives)
                .setReference(reference)
                .setRescueWindowSize(rescueWindowSize)
                .build()
                .tryRescue(position);
          }
        });
  }

  @Test
  public void testEmptyWindow() throws Exception {
    assertEquals(
        Optional.<SequenceRescuer.RescuedVariants>absent(),
        tryRescue(
            "chr1",
            10049,
            variants(variant("chr1", 8000, "G", "C", "1/1")),
            variants(variant("chr1", 10049, "CTTAAGCT", "C", "1/1")),
            variants(),
            50));
  }

  @Test
  public void testTooManyPaths() throws Exception {
    assertEquals(
        Optional.<SequenceRescuer.RescuedVariants>absent(),
        tryRescue(
            "chr1",
            10000,
            variants(
                variant("chr1", 10049, "CTTAAGCT", "C", "1/1"),
                variant("chr1", 10053, "TGCGT", "T", "0/1"),
                variant("chr1", 10055, "GCTAA", "G", "0/1"),
                variant("chr1", 10057, "TA", "T", "1/1"),
                variant("chr1", 10058, "GC", "G", "1/1")),
            variants(
                variant("chr1", 10025, "CTTAAGCT", "C", "1/1"),
                variant("chr1", 10028, "TGCGT", "T", "0/1"),
                variant("chr1", 10029, "GCTAA", "G", "0/1"),
                variant("chr1", 10032, "TA", "T", "1/1")),
            variants(),
            50));
  }

  @Test
  public void testVariantWithMismatchedRef() throws Exception {
    assertEquals(
        Optional.<SequenceRescuer.RescuedVariants>absent(),
        tryRescue(
            "chr1",
            2,
            variants(
                variant("chr2", 2, "TGC", "TAT", "1/1")),
            variants(
                variant("chr2", 3, "G", "C", "1/1"),
                variant("chr2", 4, "C", "T", "0/1")),
            variants(),
            50));
  }

  @Test
  public void testWindowTooBig() throws Exception {
    assertEquals(
        Optional.<SequenceRescuer.RescuedVariants>absent(),
        tryRescue(
            "chr1",
            10049,
            variants(
                variant("chr1", 7001, multiply(300, "ATTGTTCATGA"), "A", "1/1"),
                variant("chr1", 10100, multiply(300, "GCCTAGGGTCA"), "G", "0/1")),
            variants(variant("chr1", 10049, "CTTAAGCT", "C", "1/1")),
            variants(),
            50));
  }
}

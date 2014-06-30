package edu.berkeley.cs.amplab.smash4j;

import static edu.berkeley.cs.amplab.smash4j.SequenceRescuer.HIGHEST_END;
import static edu.berkeley.cs.amplab.smash4j.SequenceRescuer.LOWEST_START;
import static edu.berkeley.cs.amplab.smash4j.SequenceRescuer.getChoppedVariant;
import static edu.berkeley.cs.amplab.smash4j.SequenceRescuer.getSequence;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.getReference;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variant;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variants;
import static org.junit.Assert.assertEquals;

import com.google.common.base.Optional;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
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
      final NavigableMap<Integer, VariantProto> truePositives) throws Exception {
    return FastaReader.create(reference).read(
        new FastaReader.Callback<Optional<SequenceRescuer.RescuedVariants>>() {
          @Override public Optional<SequenceRescuer.RescuedVariants> read(
              Map<String, Integer> info,
              FastaReader.Callback.FastaFile reference) {
            return SequenceRescuer.builder()
                .setContig(contig)
                .setTruePositives(truePositives)
                .setFalsePositives(falsePositives)
                .setFalseNegatives(falseNegatives)
                .setReference(reference)
                .setRescueWindowSize(50)
                .setGetTypeFunction(VariantType.getType(50))
                .build()
                .tryRescue(position);
          }
        });
  }

  @Test
  public void testEmptyWindow() throws Exception {
    assertEquals(
        Optional.absent(),
        tryRescue(
            "chr1",
            10049,
            variants(variant("chr1", 8000, "G", "C", "1/1")),
            variants(variant("chr1", 10049, "CTTAAGCT", "C", "1/1")),
            variants()));
  }

  @Test
  public void testEnlargeBounds() {
    assertEquals(SequenceRescuer.Window.create(88000, 88021), SequenceRescuer.Window.Factory.Builder
        .windowEnlarger(variants(variant("chr19", 88012, "CTTAAGCT", "C", "1/1"))).apply(
        SequenceRescuer.Window.create(88000, 88020)));
  }

  @Test
  public void testFullRescue() throws Exception {
    NavigableMap<Integer, VariantProto>
        falseNegatives = variants(
            variant("chr2", 2, "TGC", "TAT", "1/1")),
        falsePositives = variants(
            variant("chr2", 3, "G", "A", "1/1"),
            variant("chr2", 4, "C", "T", "1/1")),
        truePositives = variants();
    assertEquals(
        Optional.of(SequenceRescuer.RescuedVariants.create(falseNegatives, falsePositives)),
        tryRescue("chr2", 2, falseNegatives, falsePositives, truePositives));
    falseNegatives = variants(
        variant("chr2", 3, "GCCG", "GCA", "1/1"));
    falsePositives = variants(
        variant("chr2", 3, "GC", "G", "1/1"),
        variant("chr2", 6, "G", "A", "1/1"));
    assertEquals(
        Optional.of(SequenceRescuer.RescuedVariants.create(falseNegatives, falsePositives)),
        tryRescue("chr2", 3, falseNegatives, falsePositives, truePositives));
    falseNegatives = variants(
        variant("chr4", 3, "TC", "T", "1/1"),
        variant("chr4", 8, "C", "T", "1/1"));
    falsePositives = variants(
        variant("chr4", 4, "C", "T", "1/1"),
        variant("chr4", 7, "TC", "T", "1/1"));
    truePositives = variants(
        variant("chr4", 5, "TC", "T", "1/1"));
    assertEquals(
        Optional.of(SequenceRescuer.RescuedVariants.create(falseNegatives, falsePositives)),
        tryRescue("chr4", 3, falseNegatives, falsePositives, truePositives));

  }

  @Test
  public void testGetChoppedVariant() {
    VariantProto variant1 = variant("chr19", 88012, "CTTAAGCT", "C", "1/1"), variant2 = variant1;
    NavigableMap<Integer, VariantProto> variants = variants(variant1);
    assertEquals(Optional.of(variant1), getChoppedVariant(variants, 88015, LOWEST_START));
    assertEquals(Optional.of(variant1), getChoppedVariant(variants, 88015, HIGHEST_END));
    assertEquals(Optional.absent(), getChoppedVariant(
        variants = variants(variant1 = variant("chr19", 87962, "CTTAAGCT", "C", "1/1")), 88015,
        LOWEST_START));
    assertEquals(Optional.absent(), getChoppedVariant(variants, 88015, HIGHEST_END));
    assertEquals(Optional.of(variant2), getChoppedVariant(
        variants = variants(variant("chr19", 88008, "T", "G", "1/1"), variant2), 88015,
        LOWEST_START));
    assertEquals(Optional.of(variant2), getChoppedVariant(variants, 88015, HIGHEST_END));
    assertEquals(Optional.of(variant2), getChoppedVariant(
        variants = variants(variant("chr19", 88014, "T", "G", "1/1"), variant2), 88015,
        LOWEST_START));
    assertEquals(Optional.of(variant2), getChoppedVariant(variants, 88015, HIGHEST_END));
    assertEquals(Optional.of(variant1 = variant("chr19", 88008, "ATTGCTTAACG", "A", "0/1")),
        getChoppedVariant(variants = variants(variant1, variant2), 88008, LOWEST_START));
    assertEquals(Optional.of(variant2), getChoppedVariant(variants, 88015, HIGHEST_END));
  }

  @Test
  public void testGetSequence() throws Exception {
    assertEquals(
        Optional.of("ATTCGAAAATCG"),
        FastaReader.create(reference).read(
            new FastaReader.Callback<Optional<String>>() {
              @Override public Optional<String> read(Map<String, Integer> info,
                  FastaReader.Callback.FastaFile reference) throws Exception {
                return getSequence(
                    reference,
                    "chr3",
                    SequenceRescuer.Window.create(1, 13),
                    Arrays.asList(
                        variant("chr3", 2, "TCGA", "T", "1/1"),
                        variant("chr3", 9, "A", "AAAA", "0/1")));
              }
            }));
  }

  @Test
  public void testNormalizedVariants() throws Exception {
    NavigableMap<Integer, VariantProto>
        falseNegatives = variants(
            variant("chr4", 2, "A", "ATCTC", "0/1")),
        falsePositives = variants(
            variant("chr4", 4, "C", "CTC", "0/1"),
            variant("chr4", 6, "C", "CTC", "0/1"));
    assertEquals(
        Optional.of(SequenceRescuer.RescuedVariants.create(falseNegatives, falsePositives)),
        tryRescue("chr4", 2, falseNegatives, falsePositives, variants()));
  }

  @Test
  public void testOnlySnps() throws Exception {
    assertEquals(
        Optional.absent(),
        tryRescue(
            "chr1",
            3,
            variants(
                variant("chr1", 2, "A", "C", "1/1"),
                variant("chr1", 7, "C", "T", "0/1")),
            variants(
                variant("chr1", 4, "A", "C", "1/1")),
            variants()));
  }

  @Test
  public void testOverlappingVariants() throws Exception {
    assertEquals(
        Optional.absent(),
        tryRescue(
            "chr2",
            1,
            variants(
                variant("chr2", 1, "T", "G", "1/1")),
            variants(
                variant("chr2", 7, "GA", "A", "1/1")),
            variants(
                variant("chr2", 3, "GCC", "G", "1/1"),
                variant("chr2", 4, "C", "G", "0/1"))));
  }

  @Test
  public void testTooManyPaths() throws Exception {
    assertEquals(
        Optional.absent(),
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
            variants()));
  }

  @Test
  public void testVariantWithMismatchedRef() throws Exception {
    assertEquals(
        Optional.absent(),
        tryRescue(
            "chr1",
            2,
            variants(
                variant("chr2", 2, "TGC", "TAT", "1/1")),
            variants(
                variant("chr2", 3, "G", "C", "1/1"),
                variant("chr2", 4, "C", "T", "0/1")),
            variants()));
  }

  @Test
  public void testWindowTooBig() throws Exception {
    assertEquals(
        Optional.absent(),
        tryRescue(
            "chr1",
            10049,
            variants(
                variant("chr1", 7001, multiply(300, "ATTGTTCATGA"), "A", "1/1"),
                variant("chr1", 10100, multiply(300, "GCCTAGGGTCA"), "G", "0/1")),
            variants(variant("chr1", 10049, "CTTAAGCT", "C", "1/1")),
            variants()));
  }
}

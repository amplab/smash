package edu.berkeley.cs.amplab.smash4j;

import static edu.berkeley.cs.amplab.smash4j.TestUtils.getReference;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variant;
import static edu.berkeley.cs.amplab.smash4j.TestUtils.variants;
import static edu.berkeley.cs.amplab.smash4j.VariantEvaluator.Genotype.HET;
import static edu.berkeley.cs.amplab.smash4j.VariantEvaluator.Genotype.HOM_VAR;
import static edu.berkeley.cs.amplab.smash4j.VariantEvaluator.VariantType.SNP;
import static edu.berkeley.cs.amplab.smash4j.VariantEvaluator.VariantType.SV_INSERTION;
import static org.junit.Assert.assertEquals;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.smash4j.VariantEvaluator.GenotypeConcordance;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.NavigableMap;

public class VariantEvaluatorTest {

  private static File reference;

  private static VariantEvaluator.ContigStats contigStats(String contig,
      NavigableMap<Integer, VariantProto> trueVariants,
      NavigableMap<Integer, VariantProto> predictedVariants,
      Iterable<Integer> truePositiveLocations,
      Iterable<Integer> falsePositiveLocations,
      Iterable<Integer> falseNegativeLocations,
      Iterable<Integer> incorrectPredictions,
      GenotypeConcordance concordance) {
    return VariantEvaluator.ContigStats.create(contig,
        trueVariants,
        predictedVariants,
        truePositiveLocations,
        falsePositiveLocations,
        falseNegativeLocations,
        incorrectPredictions,
        concordance,
        50);
  }

  private static Map<String, VariantEvaluator.ContigStats> evaluate(
      final Iterable<VariantProto> trueVariants,
      final Iterable<VariantProto> predictedVariants) throws Exception {
    return evaluate(trueVariants, predictedVariants, Optional.<Iterable<VariantProto>>absent());
  }

  private static Map<String, VariantEvaluator.ContigStats> evaluate(
      final Iterable<VariantProto> trueVariants,
      final Iterable<VariantProto> predictedVariants,
      final Iterable<VariantProto> knownFalsePositives) throws Exception {
    return evaluate(trueVariants, predictedVariants, Optional.of(knownFalsePositives));
  }

  private static Map<String, VariantEvaluator.ContigStats> evaluate(
      final Iterable<VariantProto> trueVariants,
      final Iterable<VariantProto> predictedVariants,
      final Optional<Iterable<VariantProto>> knownFalsePositives) throws Exception {
    return FastaReader.create(reference).read(
        new FastaReader.Callback<Map<String, VariantEvaluator.ContigStats>>() {
          @Override public Map<String, VariantEvaluator.ContigStats> read(Map<String, Integer> info,
              FastaReader.Callback.FastaFile reference) throws Exception {
            return VariantEvaluator.builder()
                .setMaxSvBreakpointDistance(100)
                .setMaxVariantLengthDifference(100)
                .setReference(reference)
                .setRescueWindowSize(50)
                .build()
                .evaluate(trueVariants, predictedVariants, knownFalsePositives);
          }
        });
  }

  @BeforeClass
  public static void setUp() throws IOException {
    reference = getReference();
  }

  @Test
  public void testChromEvaluateGenotypeConcordance() throws Exception {
    NavigableMap<Integer, VariantProto>
        trueVariants = variants(
            variant("chr1", 2, "A", "T", "0/1"),
            variant("chr1", 5, "C", "T", "0/1"),
            variant("chr1", 9, "A", "G", "1/1")),
        predictedVariants = variants(
            variant("chr1", 2, "A", "T", "1/1"),
            variant("chr1", 6, "C", "G", "0/1"),
            variant("chr1", 9, "A", "G", "1/1"));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants,
                predictedVariants,
                Arrays.asList(2, 9),
                Arrays.asList(6),
                Arrays.asList(5),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create()
                    .increment(SNP, HET, HOM_VAR)
                    .increment(SNP, HOM_VAR, HOM_VAR))),
        evaluate(trueVariants.values(), predictedVariants.values()));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants = variants(
                    variant("chr1", 2, "A", "T", "0|1"),
                    variant("chr1", 9, "A", "G", "1|1")),
                predictedVariants = variants(
                    variant("chr1", 2, "A", "T", "1|0"),
                    variant("chr1", 9, "A", "G", "1|1")),
                Arrays.asList(2, 9),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create()
                    .increment(SNP, HET, HET)
                    .increment(SNP, HOM_VAR, HOM_VAR))),
        evaluate(trueVariants.values(), predictedVariants.values()));
  }

  @Test
  public void testChromEvaluateVariantsKnownFP() throws Exception {
    NavigableMap<Integer, VariantProto>
        trueVariants = variants(
            variant("chr1", 2, "A", "T", "0/1")),
        predictedVariants = variants(
            variant("chr1", 2, "A", "T", "0/1"),
            variant("chr1", 4, "G", "C", "1/1"),
            variant("chr1", 7, "G", "A", "0/1")),
        knownFalsePositives = variants(
            variant("chr1", 1, "A", "T", "./."),
            variant("chr1", 7, "G", Collections.<String>emptyList(), "0/0"));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                    "chr1",
                    trueVariants,
                    predictedVariants,
                    Arrays.asList(2),
                    Arrays.asList(4, 7),
                    Collections.<Integer>emptyList(),
                    Collections.<Integer>emptyList(),
                    VariantEvaluator.GenotypeConcordance.create()
                        .increment(SNP, HET, HET))
                .setKnownFalsePositives(knownFalsePositives, Arrays.asList(7), Arrays.asList(7))),
        evaluate(trueVariants.values(), predictedVariants.values(), knownFalsePositives.values()));
  }

  @Test
  public void testChromEvaluateVariantsSV() throws Exception {
    NavigableMap<Integer, VariantProto>
        trueVariants = variants(variant("chr1", 6, "C", "CGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGAAAAA", "0/1")),
        predictedVariants = variants(variant("chr1", 6, "C", "CGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGAAAAA", "0/1"));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants,
                predictedVariants,
                Arrays.asList(6),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create()
                    .increment(SV_INSERTION, HET, HET))),
        evaluate(trueVariants.values(), predictedVariants.values()));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants,
                predictedVariants = variants(variant("chr1", 6, "C", "CGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGAAAAATGC", "0/1")),
                Arrays.asList(6),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create()
                    .increment(SV_INSERTION, HET, HET))),
        evaluate(trueVariants.values(), predictedVariants.values()));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants,
                predictedVariants = variants(variant("chr1", 4, "C", "CGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGAAAAA", "0/1")),
                Arrays.asList(6),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create()
                    .increment(SV_INSERTION, HET, HET))),
        evaluate(trueVariants.values(), predictedVariants.values()));
    assertEquals(
        ImmutableMap.of(
            "chr1",
            contigStats(
                "chr1",
                trueVariants,
                predictedVariants = variants(variant("chr1", 110, "C", "CGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGACGTGAGATGAAAAA", "0/1")),
                Collections.<Integer>emptyList(),
                Arrays.asList(110),
                Arrays.asList(6),
                Collections.<Integer>emptyList(),
                VariantEvaluator.GenotypeConcordance.create())),
        evaluate(trueVariants.values(), predictedVariants.values()));
  }

  @Test
  public void testRescueChromEvalVariants() throws Exception {
    VariantProto
        variant = variant("chr2", 3, "GCCG", "GCA", "1/1");
    NavigableMap<Integer, VariantProto>
        trueVariants = variants(variant),
        predictedVariants = variants(
            variant("chr2", 3, "GC", "G", "1/1"),
            variant("chr2", 6, "G", "A", "1/1")),
        rescuedVariants = variants(variant);
    assertEquals(
        ImmutableMap.of(
            "chr2",
            contigStats(
                    "chr2",
                    trueVariants,
                    predictedVariants,
                    Collections.<Integer>emptyList(),
                    Collections.<Integer>emptyList(),
                    Collections.<Integer>emptyList(),
                    Arrays.asList(3),
                    VariantEvaluator.GenotypeConcordance.create())
                .setRescuedVariants(rescuedVariants)),
        evaluate(trueVariants.values(), predictedVariants.values()));
  }
}

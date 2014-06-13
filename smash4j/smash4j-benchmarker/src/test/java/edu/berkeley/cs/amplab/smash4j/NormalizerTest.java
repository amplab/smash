package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;

import com.google.common.base.Optional;
import com.google.common.collect.Iterables;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class NormalizerTest {

  private enum NormalizerFactory {

    CLEAN_ONLY {
      @Override
      Normalizer create(FastaReader.Callback.FastaFile fastaFile) {
        return Normalizer.cleanOnly(MAX_INDEL_SIZE, fastaFile);
      }
    },

    STANDARD {
      @Override
      Normalizer create(FastaReader.Callback.FastaFile fastaFile) {
        return Normalizer.create(MAX_INDEL_SIZE, fastaFile);
      }
    };

    private static final int MAX_INDEL_SIZE = 50;

    abstract Normalizer create(FastaReader.Callback.FastaFile fastaFile);
  }

  private static final String REFERENCE_FASTA = String.format(
      ">chr1%n" +
      "AACGCCGGA%n" +
      ">chr2%n" +
      "TTGCCGGAT%n" +
      ">chr3%n" +
      "ATCGATCGATCG%n" +
      ">chr4%n" +
      "AATCTCTCGGGG");

  private static File fastaFile;

  @BeforeClass
  public static void setUp() throws IOException {
    (fastaFile = File.createTempFile(NormalizerTest.class.getSimpleName(), ".fasta"))
        .deleteOnExit();
    try (PrintStream out = new PrintStream(new FileOutputStream(fastaFile))) {
      out.print(REFERENCE_FASTA);
    }
  }

  @Test
  public void testCleanOnly() throws Exception {
    assertEquals(
        variant("chr2", 6, "G", "CG", "0/1"),
        normalize(
            NormalizerFactory.CLEAN_ONLY,
            variant("chr2", 6, "g", "cg", "0/1")).get());
  }

  @Test
  public void testRegularRecordsAreUnchanged() throws Exception {
    VariantProto variant = variant("chr1", 2, "C", "A", "0/1");
    assertEquals(variant, normalize(variant).get());
  }

  @Test
  public void testHomRefRecordsAreRemoved() throws Exception {
    assertEquals(
        Optional.<VariantProto>absent(),
        normalize(variant("chr1", 2, "C", "C", "0/0")));
  }

  @Test
  public void testSnpsAndIndelsWithoutGenotypingAreRemoved() throws Exception {
    for (VariantProto variant : Arrays.asList(
        variant("chr1", 2, "C", "A", "."),
        variant("chr1", 3, "G", "C", "0/0"),
        variant("chr1", 4, "G", "T", "0|0"))) {
      assertEquals(Optional.<VariantProto>absent(), normalize(variant));
    }
  }

  @Test
  public void testSvWithoutGenotypingIsRetained() throws Exception {
    VariantProto variant = variant("chr1", 2, "C",
        "AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG", ".");
    assertEquals(variant, normalize(variant).get());
  }

  @Test
  public void testLowerCaseRefAltGetsUpperCased() throws Exception {
    assertEquals(
        variant("chr1", 2, "C", "A", "0/1"),
        normalize(variant("chr1", 2, "c", "a", "0/1")).get());
  }

  @Test
  public void testNormalizingAnInsertion() throws Exception {
    assertEquals(
        variant("chr1", 6, "C", "CG", "0/1", 9),
        normalize(variant("chr1", 9, "a", "ga", "0/1")).get());
  }

  @Test
  public void testNormalizingADeletion() throws Exception {
    assertEquals(
        variant("chr1", 4, "GC", "G", "0/1", 5),
        normalize(variant("chr1", 5, "cc", "c", "0/1")).get());
  }

  @Test
  public void testMultipleAltAlleles() throws Exception {
    assertEquals(
        variant("chr2", 3, "G", "GC", "0/1", 6),
        normalize(variant("chr2", 6, "G", "CG", "0/1")).get());
    VariantProto variant = variant("chr2", 3, "G", Arrays.asList("CG", "C"), "0/1");
    assertEquals(variant, normalize(variant).get());
  }

  @Test
  public void testCollidingVariants() throws Exception {
    VariantProto variant = variant("chr1", 5, "A", "TGC", "1/1");
    assertEquals(Collections.singletonList(variant),
        normalize(variant, variant("chr1", 5, "A", "GGG", "1/1")));
  }

  @Test
  public void testNormalizedToCollision() throws Exception {
    assertEquals(
        Arrays.asList(
            variant("chr2", 4, "C", "T", "0/1"),
            variant("chr2", 5, "C", "CGC", "0/1", 5),
            variant("chr4", 2, "A", "AGG", "0/1"),
            variant("chr4", 3, "T", "TCT", "0/1", 6)),
        normalize(
            variant("chr2", 4, "C", "T", "0/1"),
            variant("chr2", 5, "C", "CGC", "0/1"),
            variant("chr4", 2, "A", "AGG", "0/1"),
            variant("chr4", 6, "C", "CTC", "0/1")));
    assertEquals(
        Arrays.asList(
            variant("chr4", 2, "ATC", "A", "0/1"),
            variant("chr4", 5, "TCT", "T", "0/1", 6)),
        normalize(
            variant("chr4", 2, "ATC", "A", "0/1"),
            variant("chr4", 6, "CTC", "C", "0/1")));
  }

  @Test
  public void testNormalizeTwoToCollision() throws Exception {
    assertEquals(
        Arrays.asList(
            variant("chr4", 2, "A", "ATC", "0/1", 4),
            variant("chr4", 3, "T", "TCT", "0/1", 6)),
        normalize(
            variant("chr4", 4, "C", "CTC", "0/1"),
            variant("chr4", 6, "C", "CTC", "0/1")));
  }

  @Test
  public void testNormalizeThreeCollision() throws Exception {
    assertEquals(
        Arrays.asList(
            variant("chr4", 2, "A", "T", "0/1"),
            variant("chr4", 3, "T", "TCTTT", "0/1", 1),
            variant("chr4", 4, "CTCTC", "C", "0/1", 2)),
        normalize(
            variant("chr4", 2, "A", "ATCTT", "0/1", 1),
            variant("chr4", 2, "A", "T", "0/1"),
            variant("chr4", 2, "ATCTC", "T", "0/1", 2)));
  }

  private static VariantProto variant(String chrom, int pos, String ref, String alt,
      String genotype) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.<Integer>absent());
  }

  private static VariantProto variant(
      String chrom, int pos, String ref, List<String> alts, String genotype) {
    return variant(chrom, pos, ref, alts, genotype, Optional.<Integer>absent());
  }

  private static VariantProto variant(String chrom, int pos, String ref, String alt,
      String genotype, int originalPos) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.of(originalPos));
  }

  private static VariantProto variant(String chrom, int pos, String ref, List<String> alts,
      String genotype, Optional<Integer> originalPos) {
    VariantProto.Builder builder = VariantProto.newBuilder()
        .setContig(chrom)
        .setPosition(pos)
        .setReferenceBases(ref)
        .addAllAlternateBases(alts)
        .addCall(VariantProto.Multimap.newBuilder()
            .addEntry(VariantProto.Multimap.Entry.newBuilder()
                .setKey("GT")
                .addValue(genotype)));
    if (originalPos.isPresent()) {
      builder.getInfoBuilder().addEntry(VariantProto.Multimap.Entry.newBuilder()
          .setKey(Normalizer.NORM_INFO_TAG)
          .addValue(String.valueOf(originalPos.get())));
    }
    return builder.build();
  }

  private static Optional<VariantProto> normalize(VariantProto variant) throws Exception {
    return normalize(NormalizerFactory.STANDARD, variant);
  }

  private static Optional<VariantProto> normalize(
      final NormalizerFactory factory,
      final VariantProto variant) throws Exception {
    return Optional.fromNullable(
        Iterables.getFirst(normalize(factory, variant, new VariantProto[0]), null));
  }

  private static List<VariantProto> normalize(
      final VariantProto head,
      final VariantProto... tail) throws Exception {
    return normalize(NormalizerFactory.STANDARD,  head, tail);
  }

  private static List<VariantProto> normalize(
      final NormalizerFactory factory,
      final VariantProto head,
      final VariantProto... tail) throws Exception {
    return FastaReader.create(fastaFile).read(
        new FastaReader.Callback<List<VariantProto>>() {
          @Override public List<VariantProto> read(
              Map<String, Integer> info,
              FastaReader.Callback.FastaFile fastaFile) throws Exception {
        return factory.create(fastaFile)
            .normalize(Iterables.concat(Collections.singletonList(head), Arrays.asList(tail)))
            .toList();
          }
        });
  }
}

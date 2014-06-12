package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;

import com.google.common.base.Optional;

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
        variant("chr2", 6, "G", Collections.singletonList("CG"), "0/1"),
        normalize(
            NormalizerFactory.CLEAN_ONLY,
            variant("chr2", 6, "g", Collections.singletonList("cg"), "0/1")).get());
  }

  @Test
  public void testRegularRecordsAreUnchanged() throws Exception {
    VariantProto variant = variant("chr1", 2, "C", Collections.singletonList("A"), "0/1");
    assertEquals(variant, normalize(variant).get());
  }

  @Test
  public void testHomRefRecordsAreRemoved() throws Exception {
    assertEquals(
        Optional.<VariantProto>absent(),
        normalize(variant("chr1", 2, "C", Collections.singletonList("C"), "0/0")));
  }

  @Test
  public void testSnpsAndIndelsWithoutGenotypingAreRemoved() throws Exception {
    for (VariantProto variant : Arrays.asList(
        variant("chr1", 2, "C", Collections.singletonList("A"), "."),
        variant("chr1", 3, "G", Collections.singletonList("C"), "0/0"),
        variant("chr1", 4, "G", Collections.singletonList("T"), "0|0"))) {
      assertEquals(Optional.<VariantProto>absent(), normalize(variant));
    }
  }

  @Test
  public void testSvWithoutGenotypingIsRetained() throws Exception {
    VariantProto variant = variant("chr1", 2, "C", Collections.singletonList(
        "AAAAGAAAGGCATGACCTATCCACCCATGCCACCTGGATGGACCTCACAGGCACACTGCTTCATGAGAGAG"), ".");
    assertEquals(variant, normalize(variant).get());
  }

  @Test
  public void testLowerCaseRefAltGetsUpperCased() throws Exception {
    assertEquals(
        variant("chr1", 2, "C", Collections.singletonList("A"), "0/1"),
        normalize(variant("chr1", 2, "c", Collections.singletonList("a"), "0/1")).get());
  }

  @Test
  public void testNormalizingAnInsertion() throws Exception {
    assertEquals(
        variant("chr1", 6, "C", Collections.singletonList("CG"), "0/1", 9),
        normalize(variant("chr1", 9, "a", Collections.singletonList("ga"), "0/1")).get());
  }

  @Test
  public void testNormalizingADeletion() throws Exception {
    assertEquals(
        variant("chr1", 4, "GC", Collections.singletonList("G"), "0/1", 5),
        normalize(variant("chr1", 5, "cc", Collections.singletonList("c"), "0/1")).get());
  }

  @Test
  public void testMultipleAltAlleles() throws Exception {
    assertEquals(
        variant("chr2", 3, "G", Collections.singletonList("GC"), "0/1", 6),
        normalize(variant("chr2", 6, "G", Collections.singletonList("CG"), "0/1")).get());
    VariantProto variant = variant("chr2", 3, "G", Arrays.asList("CG", "C"), "0/1");
    assertEquals(variant, normalize(variant).get());
  }

  private static VariantProto variant(
      String chrom, int pos, String ref, List<String> alts, String genotype) {
    return variant(chrom, pos, ref, alts, genotype, Optional.<Integer>absent());
  }

  private static VariantProto variant(
      String chrom, int pos, String ref, List<String> alts, String genotype, int originalPos) {
    return variant(chrom, pos, ref, alts, genotype, Optional.of(originalPos));
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

  private static Optional<VariantProto> normalize(final NormalizerFactory factory,
      final VariantProto variant) throws Exception {
    return FastaReader.create(fastaFile).read(
        new FastaReader.Callback<Optional<VariantProto>>() {
          @Override public Optional<VariantProto> read(Map<String, Integer> info,
              FastaReader.Callback.FastaFile fastaFile) throws Exception {
            return factory.create(fastaFile).normalize(Collections.singletonList(variant)).first();
          }
        });
  }
}

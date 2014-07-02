package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableListMultimap;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;

class TestUtils {

  private static final String REFERENCE_FASTA = String.format(
      ">chr1%n" +
      "AACGCCGGA%n" +
      ">chr2%n" +
      "TTGCCGGAT%n" +
      ">chr3%n" +
      "ATCGATCGATCG%n" +
      ">chr4%n" +
      "AATCTCTCGGGG");

  static File getReference() throws IOException {
    File fastaFile = File.createTempFile(NormalizerTest.class.getSimpleName(), ".fasta");
    try (PrintStream out = new PrintStream(new FileOutputStream(fastaFile))) {
      out.print(REFERENCE_FASTA);
    }
    fastaFile.deleteOnExit();
    return fastaFile;
  }

  static Variant variant(
      String chrom, int pos, String ref, List<String> alts, String genotype) {
    return variant(chrom, pos, ref, alts, genotype, Optional.<Integer>absent());
  }

  private static Variant variant(String chrom, int pos, String ref, List<String> alts,
      String genotype, Optional<Integer> originalPos) {
    Variant.Builder builder = Variant.builder()
        .setContig(chrom)
        .setPosition(pos)
        .setReferenceBases(ref)
        .setAlternateBases(alts)
        .setCalls(
            Collections.singletonList(
                Variant.Call.builder()
                    .setGenotype(Variant.GENOTYPE_SPLITTER.apply(genotype))
                    .build()));
    if (originalPos.isPresent()) {
      builder.setInfo(
          ImmutableListMultimap.<String, String>builder()
              .put(Normalizer.NORM_INFO_TAG, String.valueOf(originalPos.get()))
              .build());
    }
    return builder.build();
  }

  static Variant variant(String chrom, int pos, String ref, String alt,
      String genotype) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.<Integer>absent());
  }

  static Variant variant(String chrom, int pos, String ref, String alt,
      String genotype, int originalPos) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.of(originalPos));
  }

  static NavigableMap<Integer, Variant> variants(Variant... variants) {
    NavigableMap<Integer, Variant> result = new TreeMap<>();
    for (Variant variant : variants) {
      result.put(variant.position(), variant);
    }
    return result;
  }
}

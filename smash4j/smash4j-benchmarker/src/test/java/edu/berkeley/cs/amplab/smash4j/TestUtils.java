package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

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

  static VariantProto variant(
      String chrom, int pos, String ref, List<String> alts, String genotype) {
    return variant(chrom, pos, ref, alts, genotype, Optional.<Integer>absent());
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

  static VariantProto variant(String chrom, int pos, String ref, String alt,
      String genotype) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.<Integer>absent());
  }

  static VariantProto variant(String chrom, int pos, String ref, String alt,
      String genotype, int originalPos) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.of(originalPos));
  }

  static NavigableMap<Integer, VariantProto> variants(VariantProto... variants) {
    NavigableMap<Integer, VariantProto> result = new TreeMap<>();
    for (VariantProto variant : variants) {
      result.put(variant.getPosition(), variant);
    }
    return result;
  }
}

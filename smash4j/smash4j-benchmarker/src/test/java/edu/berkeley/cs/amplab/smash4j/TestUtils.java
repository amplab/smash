package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.Collections;
import java.util.List;

class TestUtils {

  static VariantProto variant(String chrom, int pos, String ref, String alt,
      String genotype) {
    return variant(chrom,
        pos,
        ref,
        Collections.singletonList(alt),
        genotype,
        Optional.<Integer>absent());
  }

  static VariantProto variant(
      String chrom, int pos, String ref, List<String> alts, String genotype) {
    return variant(chrom, pos, ref, alts, genotype, Optional.<Integer>absent());
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
}

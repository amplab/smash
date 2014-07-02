package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import edu.berkeley.cs.amplab.smash4j.vcf.VcfRecord;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VariantTest {

  private static final String CONTIG = "chr1";
  private static final int POSITION = 10;
  private static final String REFERENCE_BASES = "A";
  private static final List<String> ALTERNATE_BASES = Collections.singletonList("C");
  private static final List<Integer> GENOTYPE = Arrays.asList(0, 1);
  private static final Variant EXAMPLE =
      Variant.builder()
          .setContig(CONTIG)
          .setPosition(POSITION)
          .setReferenceBases(REFERENCE_BASES)
          .setAlternateBases(ALTERNATE_BASES)
          .setCalls(
              Collections.singletonList(
                  Variant.Call.builder()
                      .setGenotype(GENOTYPE)
                      .build()))
          .build();

  @Test
  public void testCreateFromVariant() {
    assertEquals(
        EXAMPLE,
        Variant.create(
            new com.google.api.services.genomics.model.Variant()
                .setContig(CONTIG)
                .setPosition((long) POSITION)
                .setReferenceBases(REFERENCE_BASES)
                .setAlternateBases(ALTERNATE_BASES)
                .setCalls(
                    Collections.singletonList(
                        new com.google.api.services.genomics.model.Call()
                            .setGenotype(GENOTYPE)))));
  }

  @Test
  public void testCreateFromVcfRecord() {
    assertEquals(
        EXAMPLE,
        Variant.create(
            VcfRecord.builder()
                .setChrom(CONTIG)
                .setPos(POSITION)
                .setRef(REFERENCE_BASES)
                .setAlt(ALTERNATE_BASES)
                .setQual(10.0)
                .setFormat(Collections.singletonList("GT"))
                .addSample(Collections.singletonList("0/1"))
                .build()));
  }
}

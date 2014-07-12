package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CommandLineTest {

  @Test
  public void testCommandLine() {
    assertEquals(
        CommandLine.builder()
            .setLhsSampleId("lhs_sample_id")
            .setLhsVcf("lhs_vcf")
            .setReferenceFai("reference_fai")
            .setReferenceFasta("reference_fasta")
            .setRhsSampleId("rhs_sample_id")
            .setRhsVcf("rhs_vcf")
            .build(),
        CommandLine.parse(
            "--lhs_sample_id=lhs_sample_id",
            "--lhs_vcf=lhs_vcf",
            "--reference_fai=reference_fai",
            "--reference_fasta=reference_fasta",
            "--rhs_sample_id=rhs_sample_id",
            "--rhs_vcf=rhs_vcf"));
  }
}

package edu.berkeley.cs.amplab.calldiff;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.CommandLine;

/**
 * Unit test for {@link CommandLine}
 */
public class CommandLineTest {

  @Test
  public void testCommandLine() {
    assertEquals(
        CommandLine.builder()
            .setLhsSampleId("lhs_sample_id")
            .setLhsVcf("lhs_vcf")
            .setPresorted(true)
            .setReferenceFai("reference_fai")
            .setReferenceFasta("reference_fasta")
            .setRhsSampleId("rhs_sample_id")
            .setRhsVcf("rhs_vcf")
            .build(),
        CommandLine.parse(
            "--lhs_sample_id=lhs_sample_id",
            "--lhs_vcf=lhs_vcf",
            "--presorted",
            "--reference_fai=reference_fai",
            "--reference_fasta=reference_fasta",
            "--rhs_sample_id=rhs_sample_id",
            "--rhs_vcf=rhs_vcf"));
  }
}

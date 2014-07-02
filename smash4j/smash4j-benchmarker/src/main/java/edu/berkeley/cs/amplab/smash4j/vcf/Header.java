package edu.berkeley.cs.amplab.smash4j.vcf;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;

import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * The header record of a VCF file. All this is essentially is a list of sample IDs.
 */
public final class Header {

  /**
   * A builder class for {@code Header} objects.
   */
  public static final class Builder {

    private final ImmutableList.Builder<String> sampleIds;

    Builder(List<String> sampleIds) {
      this.sampleIds = ImmutableList.<String>builder().addAll(sampleIds);
    }

    /**
     * Add a new sample Id to the list of sample Ids in the header.
     */
    public Builder addSampleId(String sampleId) {
      this.sampleIds.add(sampleId);
      return this;
    }

    /**
     * Build the header.
     */
    public Header build() {
      return new Header(sampleIds.build());
    }
  }

  /**
   * Static factory method for the {@code Header.Builder} class.
   */
  public static Builder builder() {
    return new Builder(Collections.<String>emptyList());
  }

  private final List<String> sampleIds;

  Header(List<String> sampleIds) {
    this.sampleIds = sampleIds;
  }

  @Override
  public boolean equals(Object obj) {
    return this == obj
        || null != obj
        && Header.class == obj.getClass()
        && Objects.equals(sampleIds(), ((Header) obj).sampleIds());
  }

  /**
   * Return the list of sample Ids.
   */
  public List<String> sampleIds() {
    return sampleIds;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(sampleIds());
  }

  /**
   * Convert this {@code Header} into its {@code Header.Builder} representation.
   */
  public Builder toBuilder() {
    return new Builder(sampleIds());
  }

  @Override
  public String toString() {
    List<String> sampleIds = sampleIds();
    return String.format(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO%s",
        sampleIds.isEmpty() ? "" : String.format("\tFORMAT\t%s", Joiner.on('\t').join(sampleIds)));
  }
}
package edu.berkeley.cs.amplab.vcfparser;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;

import java.util.List;
import java.util.Objects;

/**
 * The header record of a VCF file. All this is essentially is a list of sample IDs.
 */
public final class Header {

  private final List<String> sampleIds;

  private Header(List<String> sampleIds) {
    this.sampleIds = sampleIds;
  }

  /**
   * Static factory method for the {@code Header.Builder} class.
   */
  public static Builder builder() {
    return new Builder();
  }

  /**
   * Convert this {@code Header} into its {@code Header.Builder} representation.
   */
  public Builder toBuilder() {
    return new Builder(sampleIds());
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

  @Override
  public boolean equals(Object obj) {
    return this == obj
        || null != obj
        && Header.class == obj.getClass()
        && Objects.equals(sampleIds(), ((Header) obj).sampleIds());
  }

  @Override
  public String toString() {
    List<String> sampleIds = sampleIds();
    return String.format(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO%s",
        sampleIds.isEmpty() ? "" : String.format("\tFORMAT\t%s", Joiner.on('\t').join(sampleIds)));
  }

  /**
   * A builder class for {@code Header} objects.
   */
  public static final class Builder {

    private final ImmutableList.Builder<String> sampleIds = ImmutableList.builder();

    private Builder() {}

    private Builder(List<String> sampleIds) {
      this.sampleIds.addAll(sampleIds);
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
}

package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.model.Variant;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;

import java.util.List;
import java.util.Map;
import java.util.Objects;

public abstract class VariantAdapter {

  public static VariantAdapter fromVariant(final Variant variant) {
    return
        new VariantAdapter() {

          @Override public List<String> alternateBases() {
            return variant.getAlternateBases();
          }

          @Override public String contig() {
            return variant.getContig();
          }

          @Override public Map<String, List<String>> info() {
            return variant.getInfo();
          }

          @Override public List<String> names() {
            return variant.getNames();
          }

          @Override public long position() {
            return variant.getPosition();
          }

          @Override public String referenceBases() {
            return variant.getReferenceBases();
          }

          @Override public String toString() {
            return variant.toString();
          }
        };
  }

  public static VariantAdapter fromVcfRecord(final VcfRecord record) {
    return
        new VariantAdapter() {

          @Override public List<String> alternateBases() {
            return record.alt();
          }

          @Override public String contig() {
            return record.chrom();
          }

          @Override public Map<String, List<String>> info() {
            return record.info();
          }

          @Override public List<String> names() {
            return record.ids();
          }

          @Override public long position() {
            return record.pos();
          }

          @Override public String referenceBases() {
            return record.ref();
          }

          @Override public String toString() {
            return record.toString();
          }
        };
  }

  public abstract List<String> alternateBases();
  public abstract String contig();
  public abstract Map<String, List<String>> info();
  public abstract List<String> names();
  public abstract long position();
  public abstract String referenceBases();
  @Override public abstract String toString();

  @Override
  public final int hashCode() {
    return Objects.hash(
        alternateBases(),
        contig(),
        info(),
        names(),
        position(),
        referenceBases());
  }

  @Override
  public final boolean equals(Object obj) {
    boolean same = this == obj;
    if (!same && obj instanceof VariantAdapter) {
      VariantAdapter rhs = (VariantAdapter) obj;
      return Objects.equals(alternateBases(), rhs.alternateBases())
          && Objects.equals(contig(), rhs.contig())
          && Objects.equals(info(), rhs.info())
          && Objects.equals(names(), rhs.names())
          && Objects.equals(position(), rhs.position())
          && Objects.equals(referenceBases(), rhs.referenceBases());
    }
    return same;
  }
}

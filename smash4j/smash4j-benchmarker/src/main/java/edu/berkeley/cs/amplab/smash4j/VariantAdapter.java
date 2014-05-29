package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableMap;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public abstract class VariantAdapter {

  private static final Function<Call, Map<String, List<String>>> CALL_GET_INFO =
      new Function<Call, Map<String, List<String>>>() {
        @Override public Map<String, List<String>> apply(Call call) {
          return call.getInfo();
        }
      };

  public static VariantAdapter fromVariant(final Variant variant) {
    return
        new VariantAdapter() {

          @Override public List<String> names() {
            return variant.getNames();
          }

          @Override public String contig() {
            return variant.getContig();
          }

          @Override public long position() {
            return variant.getPosition();
          }

          @Override public String referenceBases() {
            return variant.getReferenceBases();
          }

          @Override public List<String> alternateBases() {
            return variant.getAlternateBases();
          }

          @Override public Map<String, List<String>> info() {
            return variant.getInfo();
          }

          @Override public List<Map<String, List<String>>> calls() {
            return FluentIterable.from(variant.getCalls()).transform(CALL_GET_INFO).toList();
          }

          @Override public String toString() {
            return variant.toString();
          }
        };
  }

  public static VariantAdapter fromVcfRecord(final VcfRecord record) {
    return
        new VariantAdapter() {

          private final Function<List<String>, Map<String, List<String>>> convertCallToMap =
              new Function<List<String>, Map<String, List<String>>>() {
                @Override public Map<String, List<String>> apply(List<String> call) {
                  ImmutableMap.Builder<String, List<String>> map = ImmutableMap.builder();
                  for (
                      Iterator<String> keys = record.format().iterator(), values = call.iterator();
                      keys.hasNext() && values.hasNext();
                      map.put(keys.next(), Arrays.asList(values.next().split(","))));
                  return map.build();
                }
              };

          @Override public List<String> names() {
            return record.ids();
          }

          @Override public String contig() {
            return record.chrom();
          }

          @Override public long position() {
            return record.pos();
          }

          @Override public String referenceBases() {
            return record.ref();
          }

          @Override public List<String> alternateBases() {
            return record.alt();
          }

          @Override public Map<String, List<String>> info() {
            return record.info();
          }

          @Override public List<Map<String, List<String>>> calls() {
            return FluentIterable.from(record.samples()).transform(convertCallToMap).toList();
          }

          @Override public String toString() {
            return record.toString();
          }
        };
  }

  public abstract List<String> alternateBases();

  public abstract List<Map<String, List<String>>> calls();

  public abstract String contig();

  public abstract Map<String, List<String>> info();

  public abstract List<String> names();

  public abstract long position();

  public abstract String referenceBases();

  @Override
  public final int hashCode() {
    return Objects.hash(
        names(),
        contig(),
        position(),
        referenceBases(),
        alternateBases(),
        info(),
        calls());
  }

  @Override
  public final boolean equals(Object obj) {
    boolean same = this == obj;
    if (!same && obj instanceof VariantAdapter) {
      VariantAdapter rhs = (VariantAdapter) obj;
      return Objects.equals(names(), rhs.names())
          && Objects.equals(contig(), rhs.contig())
          &&Objects.equals(position(), rhs.position())
          &&Objects.equals(referenceBases(), rhs.referenceBases())
          &&Objects.equals(alternateBases(), rhs.alternateBases())
          &&Objects.equals(info(), rhs.info())
          &&Objects.equals(calls(), rhs.calls());
    }
    return same;
  }

  @Override
  public abstract String toString();
}

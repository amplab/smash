package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;

import com.google.api.client.repackaged.com.google.common.base.Joiner;
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class VariantProtoConverterTest {

  private abstract static class ReverseVariantProtoConverter<X> {

    static Map<String, List<String>> extractMap(VariantProto.Multimap multimap) {
      ImmutableMap.Builder<String, List<String>> map = ImmutableMap.builder();
      for (VariantProto.Multimap.Entry entry : multimap.getEntryList()) {
        map.put(entry.getKey(), entry.getValueList());
      }
      return map.build();
    }

    abstract X convert(VariantProto proto);
  }

  private static final ReverseVariantProtoConverter<VcfRecord> TO_VCF_RECORD =
      new ReverseVariantProtoConverter<VcfRecord>() {
        @Override VcfRecord convert(VariantProto proto) {
          VcfRecord.Builder record = VcfRecord.builder()
              .setIds(proto.getNamesList())
              .setChrom(proto.getContig())
              .setPos(Long.valueOf(proto.getPosition()).intValue())
              .setRef(proto.getReferenceBases())
              .setAlt(proto.getAlternateBasesList())
              .setInfo(extractMap(proto.getInfo()));
          ImmutableSet.Builder<String> format = ImmutableSet.builder();
          for (VariantProto.Multimap call : proto.getCallList()) {
            ImmutableList.Builder<String> sample = ImmutableList.builder();
            for (VariantProto.Multimap.Entry entry : call.getEntryList()) {
              format.add(entry.getKey());
              sample.add(Joiner.on(",").join(entry.getValueList()));
            }
            record.addSample(sample.build());
          }
          return record.setFormat(ImmutableList.copyOf(format.build())).build();
        }
      };

  private static final ReverseVariantProtoConverter<Variant> TO_VARIANT =
      new ReverseVariantProtoConverter<Variant>() {

        private final Function<VariantProto.Multimap, Call> callFromMultimap =
            new Function<VariantProto.Multimap, Call>() {
              @Override public Call apply(VariantProto.Multimap multimap) {
                return new Call().setInfo(extractMap(multimap));
              }
            };

        @Override Variant convert(VariantProto proto) {
          return new Variant()
              .setNames(proto.getNamesList())
              .setContig(proto.getContig())
              .setPosition(Long.valueOf(proto.getPosition()))
              .setReferenceBases(proto.getReferenceBases())
              .setAlternateBases(proto.getAlternateBasesList())
              .setInfo(extractMap(proto.getInfo()))
              .setCalls(
                  FluentIterable.from(proto.getCallList()).transform(callFromMultimap).toList());
        }
      };

  @Test
  public void testVcfRecordConverter() {
    testConverter(VariantProtoConverter.VCF_RECORD_CONVERTER, TO_VCF_RECORD);
  }

  @Test
  public void testVariantConverter() {
    testConverter(VariantProtoConverter.VARIANT_CONVERTER, TO_VARIANT);
  }

  private static final VariantProto EXAMPLE_PROTO = VariantProto.newBuilder()
      .addNames("rs149201999")
      .setContig("22")
      .setPosition(16050408)
      .setReferenceBases("T")
      .addAlternateBases("C")
      .setInfo(multimap(
          entry("LDAF", "0.0649"),
          entry("RSQ", "0.8652"),
          entry("AN", "2184"),
          entry("ERATE", "0.0046"),
          entry("VT", "SNP"),
          entry("AA", "."),
          entry("AVGPOST", "0.9799"),
          entry("THETA", "0.0149"),
          entry("SNPSOURCE", "LOWCOV"),
          entry("AC", "134"),
          entry("AF", "0.06"),
          entry("ASN_AF", "0.04"),
          entry("AMR_AF", "0.05"),
          entry("AFR_AF", "0.10"),
          entry("EUR_AF", "0.06")))
      .addCall(multimap(
          entry("GT", "0|0"),
          entry("DS", "0.050"),
          entry("GL", "-0.03", "-1.17", "-5.00")))
      .build();

  private static VariantProto.Multimap multimap(VariantProto.Multimap.Entry... entries) {
    return VariantProto.Multimap.newBuilder()
        .addAllEntry(Arrays.asList(entries))
        .build();
  }

  private static VariantProto.Multimap.Entry entry(String key, String... values) {
    return VariantProto.Multimap.Entry.newBuilder()
        .setKey(key)
        .addAllValue(Arrays.asList(values))
        .build();
  }

  private static <X> void testConverter(
      VariantProtoConverter<X> converter,
      ReverseVariantProtoConverter<X> reverseConverter) {
    assertEquals(EXAMPLE_PROTO, converter.convert(reverseConverter.convert(EXAMPLE_PROTO)));
  }
}

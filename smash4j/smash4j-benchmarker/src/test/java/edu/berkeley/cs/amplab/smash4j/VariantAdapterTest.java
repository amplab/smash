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
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

@RunWith(JUnit4.class)
public class VariantAdapterTest {

  private static final Function<Map<String, List<String>>, Call> INFO_TO_CALL =
      new Function<Map<String, List<String>>, Call>() {
        @Override public Call apply(Map<String, List<String>> info) {
          return new Call().setInfo(info);
        }
      };

  @Test
  public void testVariantAdapter() {
    testEqualsAndHashCode(
        Collections.singletonList("C"),
        Collections.<Map<String, List<String>>>singletonList(
            ImmutableMap.<String, List<String>>builder()
                .put("GT", Collections.singletonList("0|0"))
                .put("DS", Collections.singletonList("0.050"))
                .put("GL", Arrays.asList("-0.03", "-1.17", "-5.00"))
                .build()),
        "22",
        ImmutableMap.<String, List<String>>builder()
            .put("LDAF", Collections.singletonList("0.0649"))
            .put("RSQ", Collections.singletonList("0.8652"))
            .put("AN", Collections.singletonList("2184"))
            .put("ERATE", Collections.singletonList("0.0046"))
            .put("VT", Collections.singletonList("SNP"))
            .put("AA", Collections.singletonList("."))
            .put("AVGPOST", Collections.singletonList("0.9799"))
            .put("THETA", Collections.singletonList("0.0149"))
            .put("SNPSOURCE", Collections.singletonList("LOWCOV"))
            .put("AC", Collections.singletonList("134"))
            .put("AF", Collections.singletonList("0.06"))
            .put("ASN_AF", Collections.singletonList("0.04"))
            .put("AMR_AF", Collections.singletonList("0.05"))
            .put("AFR_AF", Collections.singletonList("0.10"))
            .put("EUR_AF", Collections.singletonList("0.06"))
            .build(),
        Collections.singletonList("rs149201999"),
        16050408,
        "T");
  }

  private static void testEqualsAndHashCode(
      List<String> alt,
      List<Map<String, List<String>>> calls,
      String contig,
      Map<String, List<String>> info,
      List<String> names,
      int pos,
      String ref) {
    VcfRecord.Builder record = VcfRecord.builder()
        .setChrom(contig)
        .setInfo(info)
        .setIds(names)
        .setPos(pos)
        .setRef(ref);
    final ImmutableSet.Builder<String> format = ImmutableSet.builder();
    for (Map<String, List<String>> map : calls) {
      ImmutableList.Builder<String> sample = ImmutableList.builder();
      for (Map.Entry<String, List<String>> entry : map.entrySet()) {
        String key = entry.getKey();
        List<String> value = entry.getValue();
        format.add(key);
        sample.add(Joiner.on(',').join(value));
      }
      record.addSample(sample.build());
    }
    record.setFormat(ImmutableList.copyOf(format.build()));
    VariantAdapter
        fromVariant = VariantAdapter.fromVariant(new Variant()
            .setCalls(FluentIterable.from(calls).transform(INFO_TO_CALL).toList())
            .setContig(contig)
            .setInfo(info)
            .setNames(names)
            .setPosition((long) pos)
            .setReferenceBases(ref)),
        fromVcfRecord = VariantAdapter.fromVcfRecord(record.build());
    assertEquals(fromVariant, fromVcfRecord);
    assertEquals(fromVcfRecord, fromVariant);
    assertEquals(fromVariant.hashCode(), fromVcfRecord.hashCode());
  }

}

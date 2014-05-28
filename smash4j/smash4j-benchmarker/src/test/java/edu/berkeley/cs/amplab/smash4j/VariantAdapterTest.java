package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;

import com.google.api.services.genomics.model.Variant;
import com.google.common.collect.ImmutableMap;
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

  @Test
  public void testVariantAdapter() {
    testVariantAdapter(
        Collections.singletonList("A"),
        "chr20",
        ImmutableMap.of(
            "foo", Collections.singletonList("bar"),
            "fizz", Arrays.asList("buzz", "fizzbuzz")),
        Arrays.asList("hello", "goodbye"),
        123,
        "T");
  }

  private static void testVariantAdapter(
      List<String> alt,
      String contig,
      Map<String, List<String>> info,
      List<String> names,
      Integer pos,
      String ref) {
    Variant variant = new Variant();
    VcfRecord.Builder record = VcfRecord.builder();
    if (!(null == alt || alt.isEmpty())) {
      variant.setAlternateBases(alt);
      record.setAlt(alt);
    }
    if (null != contig) {
      variant.setContig(contig);
      record.setChrom(contig);
    }
    if (!(null == info || info.isEmpty())) {
      variant.setInfo(info);
      record.setInfo(info);
    }
    if (!(null == names || names.isEmpty())) {
      variant.setNames(names);
      record.setIds(names);
    }
    if (null != pos) {
      variant.setPosition(pos.longValue());
      record.setPos(pos);
    }
    if (null != ref) {
      variant.setReferenceBases(ref);
      record.setRef(ref);
    }
    VariantAdapter
        fromVariant = VariantAdapter.fromVariant(variant),
        fromVcfRecord = VariantAdapter.fromVcfRecord(record.build());
    assertEquals(fromVariant, fromVcfRecord);
    assertEquals(fromVcfRecord, fromVariant);
  }

}

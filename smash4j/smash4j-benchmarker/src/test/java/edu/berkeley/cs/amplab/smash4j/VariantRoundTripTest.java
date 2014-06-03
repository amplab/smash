package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.google.common.collect.FluentIterable;

import org.junit.Test;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class VariantRoundTripTest {

  @Test
  public void testRoundTrip() throws Exception {
    final List<VariantProto> variants = Arrays.asList(
        VariantProto.newBuilder()
            .addNames("rs6054257")
            .setContig("20")
            .setPosition(14370)
            .setReferenceBases("G")
            .addAlternateBases("A")
            .setInfo(multimap(entry("NS", "3"), entry("DP", "14"), entry("AF", "0.5"), entry("DB"),
                entry("H2")))
            .addCall(multimap(entry("GT", "0|0"), entry("GQ", "48"), entry("DP", "1"),
                entry("HQ", "51", "51")))
            .addCall(multimap(entry("GT", "1|0"), entry("GQ", "48"), entry("DP", "8"),
                entry("HQ", "51", "51")))
            .addCall(multimap(entry("GT", "1/1"), entry("GQ", "43"), entry("DP", "5"),
                entry("HQ", ".", ".")))
            .build(),
        VariantProto.newBuilder()
            .setContig("20")
            .setPosition(17330)
            .setReferenceBases("T")
            .addAlternateBases("A")
            .setInfo(multimap(entry("NS", "3"), entry("DP", "11"), entry("AF", "0.017")))
            .addCall(multimap(entry("GT", "0|1"), entry("GQ", "3"), entry("DP", "5"),
                entry("HQ 0|0", "65", "3")))
            .addCall(multimap(entry("GT", "0/0"), entry("GQ", "41"), entry("DP", "3")))
            .build(),
        VariantProto.newBuilder()
            .addNames("rs6040355")
            .setContig("20")
            .setPosition(1110696)
            .setReferenceBases("A")
            .addAlternateBases("G")
            .addAlternateBases("T")
            .setInfo(multimap(entry("NS", "2"), entry("DP", "10"), entry("AF", "0.333", "0.667"),
                entry("AA", "T"), entry("DB")))
            .addCall(multimap(entry("GT", "1|2"), entry("GQ", "21"), entry("DP", "6"),
                entry("HQ", "23", "27")))
            .addCall(multimap(entry("GT", "2|1"), entry("GQ", "2"), entry("DP", "0"),
                entry("HQ", "18", "2")))
            .addCall(multimap(entry("GT", "2/2"), entry("GQ", "35"), entry("DP", "4")))
            .build(),
        VariantProto.newBuilder()
            .setContig("20")
            .setPosition(1230237)
            .setReferenceBases("T")
            .setInfo(multimap(entry("NS", "3"), entry("DP", "13"), entry("AA", "T")))
            .addCall(multimap(entry("GT", "0|0"), entry("GQ", "48"), entry("DP", "4"),
                entry("HQ 0|0", "51", "51")))
            .addCall(multimap(entry("GT", "0/0"), entry("GQ", "61"), entry("DP", "2")))
            .build(),
        VariantProto.newBuilder()
            .addNames("microsat1")
            .setContig("20")
            .setPosition(1234567)
            .setReferenceBases("GTC")
            .addAlternateBases("G")
            .addAlternateBases("GTCT")
            .setInfo(multimap(entry("NS", "3"), entry("DP", "9"), entry("AA", "G")))
            .addCall(multimap(entry("GT", "0/1"), entry("GQ", "35"), entry("DP", "4")))
            .addCall(multimap(entry("GT", "0/2"), entry("GQ", "17"), entry("DP", "2")))
            .addCall(multimap(entry("GT", "1/1"), entry("GQ", "40"), entry("DP", "3")))
            .build());
    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
    VariantWriter.create(buffer).writeVariants(VariantScanner.fromVariantProtos(variants));
    assertEquals(
        variants,
        VariantScanner.fromVariantProtoInputStream(new ByteArrayInputStream(buffer.toByteArray()))
            .scan(
                new VariantScanner.Callback<List<VariantProto>>() {
                  @Override public List<VariantProto> scan(
                      final FluentIterable<VariantProto> actual) throws Exception {
                    return VariantScanner.fromVariantProtos(variants).scan(
                        new VariantScanner.Callback<List<VariantProto>>() {
                          @Override
                          public List<VariantProto> scan(FluentIterable<VariantProto> expected) {
                            Iterator<VariantProto>
                                lhs = expected.iterator(), rhs = actual.iterator();
                            while (lhs.hasNext() && rhs.hasNext()) {
                              assertEquals(lhs.next(), rhs.next());
                            }
                            assertFalse(lhs.hasNext() || rhs.hasNext());
                            return variants;
                          }
                        });
                  }
                }));
  }

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
}

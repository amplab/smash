package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;

public class VcfCallScannerTest {

  private static File makeTempFile(String resource, String prefix, String suffix)
      throws IOException {
    try (InputStream in = VcfCallScannerTest.class.getResourceAsStream(resource)) {
      File file = File.createTempFile(prefix, suffix);
      file.deleteOnExit();
      try (OutputStream out = new FileOutputStream(file)) {
        byte[] buffer = new byte[4096];
        for (int n = in.read(buffer); 0 < n; n = in.read(buffer)) {
          out.write(buffer, 0, n);
        }
      }
      return file;
    }
  }

  private static void testScan(VcfCallScanner scanner, Call... expected) throws IOException {
    assertEquals(
        Arrays.asList(expected),
        scanner.scan(stream -> stream.collect(Collectors.toList())));
  }

  private VcfCallScanner na00001;
  private VcfCallScanner na00002;
  private VcfCallScanner na00003;

  @Before
  public void setUp() throws IOException {
    File vcf = makeTempFile("/edu/berkeley/cs/amplab/vardiff/sample.vcf", "temp", ".vcf");
    na00001 = VcfCallScanner.create(vcf, "NA00001");
    na00002 = VcfCallScanner.create(vcf, "NA00002");
    na00003 = VcfCallScanner.create(vcf, "NA00003");
  }

  @Test
  public void testScan() throws IOException {
    testScan(
        na00001,
        TestCall.create("20", 14370, "G", Collections.singletonList("A"), Arrays.asList(0, 0),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 17330, "T", Collections.singletonList("A"), Arrays.asList(0, 0),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1110696, "A", Arrays.asList("G", "T"), Arrays.asList(1, 2),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1230237, "T", Collections.emptyList(), Arrays.asList(0, 0),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1234567, "GTC", Arrays.asList("G", "GTCT"), Arrays.asList(0, 1)));
    testScan(
        na00002,
        TestCall.create("20", 14370, "G", Collections.singletonList("A"), Arrays.asList(1, 0),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 17330, "T", Collections.singletonList("A"), Arrays.asList(0, 1),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1110696, "A", Arrays.asList("G", "T"), Arrays.asList(2, 1),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1230237, "T", Collections.emptyList(), Arrays.asList(0, 0),
            Call.Phaseset.DEFAULT),
        TestCall.create("20", 1234567, "GTC", Arrays.asList("G", "GTCT"), Arrays.asList(0, 2)));
    testScan(
        na00003,
        TestCall.create("20", 14370, "G", Collections.singletonList("A"), Arrays.asList(1, 1)),
        TestCall.create("20", 17330, "T", Collections.singletonList("A"), Arrays.asList(0, 0)),
        TestCall.create("20", 1110696, "A", Arrays.asList("G", "T"), Arrays.asList(2, 2)),
        TestCall.create("20", 1230237, "T", Collections.emptyList(), Arrays.asList(0, 0)),
        TestCall.create("20", 1234567, "GTC", Arrays.asList("G", "GTCT"), Arrays.asList(1, 1)));
  }
}

package edu.berkeley.cs.amplab.vcfparser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;

import com.google.common.collect.FluentIterable;
import org.junit.Test;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.net.URL;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class VcfParserTest {

  @Test
  public void testAllMetaInfoRecordTypes() throws IOException {
    serializeAndParseTest(
        MetaInformation.builder(MetaInformation.FileFormat.Format.V4_1)
            .addAssembly(MetaInformation.Assembly.create(new URL("http://domain.com/assembly")))
            .addAlt(MetaInformation.Alt.builder()
                .setId(MetaInformation.Alt.Type.DUP)
                .setDescription("duplication")
                .addExtraField("Ef1", "value1"))
            .addContig(MetaInformation.Contig.builder()
                .setId("contigID")
                .addExtraField("URL", "http://domain.com/contig"))
            .addFilter(MetaInformation.Filter.builder()
                .setId("filterID")
                .setDescription("Des\"rip\\ion")
                .addExtraField("Ef2", "value2"))
            .addFormat(MetaInformation.Format.builder()
                .setId("formatID")
                .setNumber(MetaInformation.Number.A)
                .setType(MetaInformation.Format.Type.INTEGER)
                .setDescription("description")
                .addExtraField("Ef3", "value3"))
            .addInfo(MetaInformation.Info.builder()
                .setId("infoID")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Info.Type.STRING)
                .setDescription("desc")
                .addExtraField("Ef5", "value5"))
            .addPedigree(MetaInformation.Pedigree.builder()
                .addExtraField("foo", "bar"))
            .addPedigreeDB(MetaInformation.PedigreeDB.create(new URL("http://domain.com/pedigreeDB")))
            .addSample(MetaInformation.Sample.builder()
                .setId("sampleID")
                .setGenome("G")
                .setMixture("M")
                .setDescription("D"))
            .addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create("foo", "bar"))
            .build(),
        Header.builder()
            .build());
  }

  @Test
  public void testExampleInVcfSpec() throws IOException {
    serializeAndParseTest(
        MetaInformation.builder(MetaInformation.FileFormat.Format.V4_2)
            .addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create("fileDate", "20090805"))
            .addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create("source", "myImputationProgramV3.1"))
            .addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create("reference", "file:///seq/references/1000GenomesPilot-NCBI36.fasta"))
            .addContig(MetaInformation.Contig.builder()
                .setId("20")
                .addExtraField("length", "62435964")
                .addExtraField("assembly", "B36")
                .addExtraField("md5", "f126cdf8a6e0c7f379d618ff66beb2da")
                .addExtraField("species", "\"Homo sapiens\"")
                .addExtraField("taxonomy", "x"))
            .addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create("phasing", "partial"))
            .addInfo(MetaInformation.Info.builder()
                .setId("NS")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Info.Type.INTEGER)
                .setDescription("Number of Samples With Data"))
            .addInfo(MetaInformation.Info.builder()
                .setId("DP")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Info.Type.INTEGER)
                .setDescription("Total Depth"))
            .addInfo(MetaInformation.Info.builder()
                .setId("AF")
                .setNumber(MetaInformation.Number.A)
                .setType(MetaInformation.Info.Type.FLOAT)
                .setDescription("Allele Frequency"))
            .addInfo(MetaInformation.Info.builder()
                .setId("AA")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Info.Type.STRING)
                .setDescription("Ancestral Allele"))
            .addInfo(MetaInformation.Info.builder()
                .setId("DB")
                .setNumber(MetaInformation.Number.create(0))
                .setType(MetaInformation.Info.Type.FLAG)
                .setDescription("dbSNP membership, build 129"))
            .addInfo(MetaInformation.Info.builder()
                .setId("H2")
                .setNumber(MetaInformation.Number.create(0))
                .setType(MetaInformation.Info.Type.FLAG)
                .setDescription("HapMap2 membership"))
            .addFilter(MetaInformation.Filter.builder()
                .setId("q10")
                .setDescription("Quality below 10"))
            .addFilter(MetaInformation.Filter.builder()
                .setId("s50")
                .setDescription("Less than 50% of samples have data"))
            .addFormat(MetaInformation.Format.builder()
                .setId("GT")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Format.Type.STRING)
                .setDescription("Genotype"))
            .addFormat(MetaInformation.Format.builder()
                .setId("GQ")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Format.Type.INTEGER)
                .setDescription("Genotype Quality"))
            .addFormat(MetaInformation.Format.builder()
                .setId("DP")
                .setNumber(MetaInformation.Number.create(1))
                .setType(MetaInformation.Format.Type.INTEGER)
                .setDescription("Read Depth"))
            .addFormat(MetaInformation.Format.builder()
                .setId("HQ")
                .setNumber(MetaInformation.Number.create(2))
                .setType(MetaInformation.Format.Type.INTEGER)
                .setDescription("Haplotype Quality"))
            .build(),
        Header.builder()
            .addSampleId("NA00001")
            .addSampleId("NA00002")
            .addSampleId("NA00003")
            .build(),
        VcfRecord.builder()
            .setChrom("20")
            .setPos(14370)
            .setId("rs6054257")
            .setRef("G")
            .setAlt("A")
            .setQual(29)
            .setFilter("PASS")
            .setInfo("NS=3;DP=14;AF=0.5;DB;H2")
            .setFormat("GT:GQ:DP:HQ")
            .addSample("0|0:48:1:51,51")
            .addSample("1|0:48:8:51,51")
            .addSample("1/1:43:5:.,.")
            .build(),
        VcfRecord.builder()
            .setChrom("20")
            .setPos(17330)
            .setRef("T")
            .setAlt("A")
            .setQual(3)
            .setFilter("q10")
            .setInfo("NS=3;DP=11;AF=0.017")
            .setFormat("GT:GQ:DP:HQ")
            .addSample("0|0:49:3:58,50")
            .addSample("0|1:3:5:65,3")
            .addSample("0/0:41:3")
            .build(),
        VcfRecord.builder()
            .setChrom("20")
            .setPos(1110696)
            .setId("rs6040355")
            .setRef("A")
            .setAlt("G,T")
            .setQual(67)
            .setFilter("PASS")
            .setInfo("NS=2;DP=10;AF=0.333,0.667;AA=T;DB")
            .setFormat("GT:GQ:DP:HQ")
            .addSample("1|2:21:6:23,27")
            .addSample("2|1:2:0:18,2")
            .addSample("2/2:35:4")
            .build(),
        VcfRecord.builder()
            .setChrom("20")
            .setPos(1230237)
            .setRef("T")
            .setQual(47)
            .setFilter("PASS")
            .setInfo("NS=3;DP=13;AA=T")
            .setFormat("GT:GQ:DP:HQ")
            .addSample("0|0:54:7:56,60")
            .addSample("0|0:48:4:51,51")
            .addSample("0/0:61:2")
            .build(),
        VcfRecord.builder()
            .setChrom("20")
            .setPos(1234567)
            .setId("microsat1")
            .setRef("GTC")
            .setAlt("G,GTCT")
            .setQual(50)
            .setFilter("PASS")
            .setInfo("NS=3;DP=9;AA=G")
            .setFormat("GT:GQ:DP")
            .addSample("0/1:35:4")
            .addSample("0/2:17:2")
            .addSample("1/1:40:3")
            .build()
    );
  }

  /**
   * This test serializes its inputs into a buffer, then parses the resulting buffer
   * back into objects, verifying that the objects that were serialized are equal to
   * the objects reified from the parse.
   */
  private static void serializeAndParseTest(
      final MetaInformation metainfo,
      final Header header,
      VcfRecord... vcfRecords) {
    final Writer buffer = new StringWriter();
    final List<VcfRecord> records = Arrays.asList(vcfRecords);
    assertSame(buffer, VcfWriter.to(buffer).write(metainfo, header, records));
    assertSame(buffer, VcfReader.from(new StringReader(buffer.toString())).read(
        new VcfReader.Callback<Writer>() {
          @Override
          public Writer readVcf(
              MetaInformation actualMetainfo,
              Header actualHeader,
              FluentIterable<VcfRecord> actualRecords) {
            assertEquals(metainfo, actualMetainfo);
            assertEquals(header, actualHeader);
            Iterator<?> lhs = records.iterator(), rhs = actualRecords.iterator();
            for (; lhs.hasNext() && rhs.hasNext(); assertEquals(lhs.next(), rhs.next()));
            assertFalse(lhs.hasNext() || rhs.hasNext());
            return buffer;
          }
        }));
  }
}

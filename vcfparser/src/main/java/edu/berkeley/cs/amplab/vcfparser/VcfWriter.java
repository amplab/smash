package edu.berkeley.cs.amplab.vcfparser;

import java.io.Closeable;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

/**
 * A writer for VCF data.
 */
public class VcfWriter implements Closeable {

  /**
   * Static factory method that writes VCF to the given {@code OutputStream}.
   */
  public static VcfWriter from(OutputStream out) {
    return from(new OutputStreamWriter(out));
  }

  /**
   * Static factory method that writes VCF to the given {@code Writer}.
   */
  public static VcfWriter from(Writer out) {
    return new VcfWriter(out instanceof PrintWriter ? (PrintWriter) out : new PrintWriter(out));
  }

  private final PrintWriter out;

  private VcfWriter(PrintWriter out) {
    this.out = out;
  }

  /**
   * Write the given {@code MetaInformation}, {@code Header}, and {@code Iterable<VcfRecord>}s to
   * the output sink.
   */
  public VcfWriter write(MetaInformation metaInformation, Header header, Iterable<VcfRecord> records) {
    out.print(metaInformation);
    out.println(header);
    for (VcfRecord record : records) {
      out.println(record);
    }
    return this;
  }

  @Override
  public void close() {
    out.close();
  }
}

package edu.berkeley.cs.amplab.smash4j.vcf;

import com.google.common.base.Function;
import com.google.common.base.Functions;

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

/**
 * A writer for VCF data.
 *
 * @param <S> The type of object (sink) being written to. Returned by the {@link #write} method.
 */
public class VcfWriter<S> {

  private static final Function<OutputStream, Writer> NEW_OUTPUT_STREAM_WRITER =
      new Function<OutputStream, Writer>() {
        @Override public Writer apply(OutputStream out) {
          return new OutputStreamWriter(out);
        }
      };

  /**
   * Static factory method that writes VCF to the given {@code OutputStream}.
   */
  public static <OS extends OutputStream> VcfWriter<OS> to(OS outputStream) {
    return create(outputStream, NEW_OUTPUT_STREAM_WRITER);
  }

  private static <S> VcfWriter<S> create(S sink, Function<? super S, Writer> writerFactory) {
    Writer out = writerFactory.apply(sink);
    return new VcfWriter<>(sink, out instanceof PrintWriter ? (PrintWriter) out : new PrintWriter(out));
  }

  /**
   * Static factory method that writes VCF to the given {@code Writer}.
   */
  public static <W extends Writer> VcfWriter<W> to(W writer) {
    return create(writer, Functions.<Writer>identity());
  }
  private final PrintWriter out;
  private final S sink;

  private VcfWriter(S sink, PrintWriter out) {
    this.sink = sink;
    this.out = out;
  }

  /**
   * Write the given metainfo, header, and records to the sink, then return the sink.
   */
  public S write(MetaInformation metaInformation, Header header, Iterable<VcfRecord> records) {
    out.print(metaInformation);
    out.println(header);
    for (VcfRecord record : records) {
      out.println(record);
    }
    out.flush();
    return sink;
  }
}

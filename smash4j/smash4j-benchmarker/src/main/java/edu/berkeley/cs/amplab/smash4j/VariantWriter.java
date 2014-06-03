package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.OutputStream;

public class VariantWriter {

  public static VariantWriter create(OutputStream out) {
    return new VariantWriter(out, Functions.<Iterable<VariantProto>>identity());
  }

  private final OutputStream out;
  private final
      Function<? super Iterable<VariantProto>, ? extends Iterable<VariantProto>> transformer;

  private VariantWriter(
      OutputStream out,
      Function<? super Iterable<VariantProto>, ? extends Iterable<VariantProto>> transformer) {
    this.out = out;
    this.transformer = transformer;
  }

  public VariantWriter transform(
      Function<? super Iterable<VariantProto>, ? extends Iterable<VariantProto>> transformer) {
    return new VariantWriter(out, Functions.compose(transformer, this.transformer));
  }

  public VariantWriter writeVariants(VariantScanner scanner) throws Exception {
    return scanner.scan(
        new VariantScanner.Callback<VariantWriter>() {
          @Override public VariantWriter scan(FluentIterable<VariantProto> variants)
              throws Exception {
            for (VariantProto variant : variants) {
              variant.writeDelimitedTo(out);
            }
            return VariantWriter.this;
          }
        });
  }
}

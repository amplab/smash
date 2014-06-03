package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.Variant;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.vcfparser.Header;
import edu.berkeley.cs.amplab.vcfparser.MetaInformation;
import edu.berkeley.cs.amplab.vcfparser.VcfReader;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;

public abstract class VariantScanner {

  public interface Callback<X> {

    X accept(FluentIterable<VariantProto> protos) throws Exception;
  }

  public static VariantScanner create(final String spec) {
    return spec.startsWith("callsets:")
        ? new VariantScanner() {
            @Override public <X> X scan(final Callback<X> callback) throws Exception {
              return VariantFetcher.create().fetchVariants(
                  Arrays.asList(spec.substring(9).split(",")),
                  new VariantFetcher.Callback<X>() {
                    @Override public X accept(Map<String, Callset> callsets,
                        FluentIterable<Variant> variants) throws Exception {
                      return callback.accept(
                          variants.transform(VariantProtoConverter.VARIANT_CONVERTER));
                    }
                  });
            }
          }
        : spec.startsWith("protofile:")
            ? new VariantScanner() {
              @Override public <X> X scan(final Callback<X> callback) throws Exception {
                  try (final InputStream in = new FileInputStream(new File(spec.substring(10)))) {
                    class ExceptionWrapper extends RuntimeException {

                      ExceptionWrapper(IOException cause) {
                        super(cause);
                      }

                      @Override public IOException getCause() {
                        return (IOException) super.getCause();
                      }
                    }
                    try {
                      return callback.accept(
                          new FluentIterable<VariantProto>() {
                            @Override public Iterator<VariantProto> iterator() {
                              return new AbstractIterator<VariantProto>() {
                                    @Override protected VariantProto computeNext() {
                                      try {
                                        return VariantProto.parseDelimitedFrom(in);
                                      } catch (IOException e) {
                                        throw new ExceptionWrapper(e);
                                      }
                                    }
                                  };
                            }
                          });
                    } catch (ExceptionWrapper e) {
                      throw e.getCause();
                    }
                  }
              }
            }
          : new VariantScanner() {
              @Override public <X> X scan(final Callback<X> callback) throws Exception {
                try (Reader in = new FileReader(new File(spec))) {
                  return VcfReader.from(in).read(new VcfReader.Callback<X>() {
                    @Override public X readVcf(MetaInformation metaInformation,
                        Header header, FluentIterable<VcfRecord> records) throws Exception {
                      return callback.accept(
                          records.transform(VariantProtoConverter.VCF_RECORD_CONVERTER));
                    }
                  });
                }
              }
            };
  }

  public abstract <X> X scan(Callback<X> callback) throws Exception;
}

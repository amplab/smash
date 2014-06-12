package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.repackaged.com.google.common.base.Throwables;
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
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public abstract class VariantScanner {

  public interface Callback<X> {
    X scan(FluentIterable<VariantProto> protos) throws Exception;
  }

  public static VariantScanner fromCallsets(final List<String> callsetIds) {
    try {
      return fromVariantFetcher(
          VariantFetcher.create(GenomicsFactory.getDefault().fromPreferences()), callsetIds);
    } catch (Exception e) {
      throw Throwables.propagate(e);
    }
  }

  public static VariantScanner fromSpec(String spec) {
    if (spec.startsWith("vcf:")) {
      return fromVcfFile(new File(spec.substring(4)));
    } else if (spec.startsWith("callsets:")) {
      return fromCallsets(Arrays.asList(spec.substring(9).split(",")));
    } else if (spec.startsWith("protofile:")) {
      return fromVariantProtoFile(new File(spec.substring(10)));
    } else {
      throw new IllegalArgumentException(String.format("Failed to parse spec: \"%s\"", spec));
    }
  }

  public static VariantScanner fromVariantFetcher(
      final VariantFetcher fetcher, final List<String> callsetIds) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            return fetcher.fetchVariants(callsetIds,
                new VariantFetcher.Callback<X>() {
                  @Override public X accept(Map<String, Callset> callsets,
                      FluentIterable<Variant> variants) throws Exception {
                    return callback.scan(
                        variants.transform(VariantProtoConverter.VARIANT_CONVERTER));
                  }
                });
          }
        };
  }

  public static VariantScanner fromVariantProtoFile(final File file) {
    return new VariantScanner() {
          @Override public <X> X scan(Callback<X> callback) throws Exception {
            try (InputStream in = new FileInputStream(file)) {
              return protoScan(in, callback);
            }
          }
        };
  }

  public static VariantScanner fromVariantProtoInputStream(final InputStream in) {
    return new VariantScanner() {
          @Override public <X> X scan(Callback<X> callback) throws Exception {
            return protoScan(in, callback);
          }
        };
  }

  public static VariantScanner fromVariantProtos(final Collection<VariantProto> variants) {
    return new VariantScanner() {
          @Override public <X> X scan(Callback<X> callback) throws Exception {
            return callback.scan(FluentIterable.from(variants));
          }
        };
  }

  public static VariantScanner fromVcfFile(final File file) {
    return new VariantScanner() {
          @Override public <X> X scan(Callback<X> callback) throws Exception {
            try (Reader in = new FileReader(file)) {
              return vcfScan(VcfReader.from(in), callback);
            }
          }
        };
  }

  public static VariantScanner fromVcfReader(final VcfReader reader) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            return vcfScan(reader, callback);
          }
        };
  }

  private static <X> X protoScan(final InputStream in, Callback<X> callback) throws Exception,
      IOException {
    class ExceptionWrapper extends RuntimeException {

      ExceptionWrapper(IOException cause) {
        super(cause);
      }

      @Override public IOException getCause() {
        return (IOException) super.getCause();
      }
    }
    try {
      return callback.scan(
          new FluentIterable<VariantProto>() {
            @Override public Iterator<VariantProto> iterator() {
              return new AbstractIterator<VariantProto>() {
                    @Override protected VariantProto computeNext() {
                      try {
                        VariantProto variant = VariantProto.parseDelimitedFrom(in);
                        return null == variant ? endOfData() : variant;
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

  private static <X> X vcfScan(final VcfReader reader, final Callback<X> callback)
      throws Exception {
    return reader.read(
        new VcfReader.Callback<X>() {
          @Override public X readVcf(MetaInformation metaInformation, Header header,
              FluentIterable<VcfRecord> records) throws Exception {
            return callback.scan(
                records.transform(VariantProtoConverter.VCF_RECORD_CONVERTER));
          }
        });
  }

  public abstract <X> X scan(Callback<X> callback) throws Exception;
}

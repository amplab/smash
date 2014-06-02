package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.Variant;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.vcfparser.Header;
import edu.berkeley.cs.amplab.vcfparser.MetaInformation;
import edu.berkeley.cs.amplab.vcfparser.VcfReader;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;

import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.util.List;
import java.util.Map;

public abstract class VariantScanner {

  public interface Callback<X> {

    X accept(FluentIterable<VariantProto> protos) throws Exception;
  }

  public static VariantScanner from(final File file) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            try (Reader in = new FileReader(file)) {
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

  public static VariantScanner from(final Genomics genomics, final List<String> callsets) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            return VariantFetcher.create(genomics).fetchVariants(callsets,
                new VariantFetcher.Callback<X>() {
                  @Override public X accept(Map<String, Callset> callsets,
                      FluentIterable<Variant> variants) throws Exception {
                    return callback.accept(
                        variants.transform(VariantProtoConverter.VARIANT_CONVERTER));
                  }
                });
          }
        };
  }

  public abstract <X> X scan(Callback<X> callback) throws Exception;
}

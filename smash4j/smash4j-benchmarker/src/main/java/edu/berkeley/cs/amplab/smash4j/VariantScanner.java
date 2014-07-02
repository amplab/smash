package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.model.Callset;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.smash4j.vcf.Header;
import edu.berkeley.cs.amplab.smash4j.vcf.MetaInformation;
import edu.berkeley.cs.amplab.smash4j.vcf.VcfReader;
import edu.berkeley.cs.amplab.smash4j.vcf.VcfRecord;

import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.util.Collections;
import java.util.Map;

public abstract class VariantScanner {

  public interface Callback<X> {
    X scan(FluentIterable<Variant> protos) throws Exception;
  }

  public static VariantScanner fromCallset(final String callsetId) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            return VariantFetcher
                .create(GenomicsFactory.getDefault().fromPreferences())
                .fetchVariants(Collections.singletonList(callsetId),
                    new VariantFetcher.Callback<X>() {
                      @Override public X accept(Map<String, Callset> callsets,
                          FluentIterable<com.google.api.services.genomics.model.Variant> variants)
                          throws Exception {
                        return callback.scan(variants.transform(Variant.CREATE_FROM_VARIANT));
                      }
                    });
          }
        };
  }

  public static VariantScanner fromVcfFile(final File file) {
    return new VariantScanner() {
          @Override public <X> X scan(final Callback<X> callback) throws Exception {
            try (Reader in = new FileReader(file)) {
              return VcfReader.from(in).read(
                  new VcfReader.Callback<X>() {
                    @Override public X readVcf(MetaInformation metaInformation, Header header,
                        FluentIterable<VcfRecord> records) throws Exception {
                      return callback.scan(records.transform(Variant.CREATE_FROM_VCF_RECORD));
                    }
                  });
            }
          }
        };
  }

  public abstract <X> X scan(Callback<X> callback) throws Exception;
}

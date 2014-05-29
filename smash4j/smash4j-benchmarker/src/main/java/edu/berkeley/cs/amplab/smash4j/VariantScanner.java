package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.Variant;
import com.google.common.collect.FluentIterable;
import edu.berkeley.cs.amplab.vcfparser.Header;
import edu.berkeley.cs.amplab.vcfparser.MetaInformation;
import edu.berkeley.cs.amplab.vcfparser.VcfReader;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.List;

public abstract class VariantScanner {

  public interface Callback<X> {
    X accumulate(X accumulator, VariantAdapter variant);
  }

  public static VariantScanner fromFile(final File file) {
    return
        new VariantScanner() {
          @Override public <X> X scan(final X initial, final Callback<X> callback)
              throws IOException {
            try (final Reader in = new FileReader(file)) {
              return VcfReader.from(in).read(
                  new VcfReader.Callback<X>() {
                    @Override public X readVcf(MetaInformation metaInformation,
                        Header header, FluentIterable<VcfRecord> records) {
                      X accumulator = initial;
                      for (VcfRecord record : records) {
                        accumulator =
                            callback.accumulate(accumulator, VariantAdapter.fromVcfRecord(record));
                      }
                      return accumulator;
                    }
                  });
            }
          }
        };
  }

  public static VariantScanner fromCallsets(
      final Genomics genomics, final List<String> callsetIds, final List<String> contigs) {
    return
        new VariantScanner() {

          private final VariantFetcher fetcher = VariantFetcher.create(genomics);

          @Override public <X> X scan(X initial, final Callback<X> callback) throws IOException {
            return fetcher.fetchVariants(
                new SearchVariantsRequest()
                    // .setDatasetId(fetcher.getDatasetId(callsetIds))
                    // TODO: unhardcode from this test dataset id.
                    .setDatasetId("376902546192")
                    .setCallsetIds(callsetIds)
                    .setStartPosition(1L)
                    .setEndPosition(Long.MAX_VALUE),
                contigs,
                initial,
                new VariantFetcher.Callback<X>() {
                  @Override public X accumulate(X accumulator, Variant variant) throws IOException {
                    return callback.accumulate(accumulator, VariantAdapter.fromVariant(variant));
                  }
                });
          }
        };
  }

  public abstract <X> X scan(X initial, Callback<X> callback) throws IOException;
}

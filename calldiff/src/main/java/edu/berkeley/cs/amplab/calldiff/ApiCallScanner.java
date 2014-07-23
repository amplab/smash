package edu.berkeley.cs.amplab.calldiff;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.common.collect.Iterables;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

/**
 * An implementation of {@link CallScanner} that retrieves calls from the cloud rather than from
 * a VCF file.
 */
public class ApiCallScanner implements CallScanner {

  public static ApiCallScanner create(Genomics genomics, String callsetId) {
    return new ApiCallScanner(genomics, callsetId);
  }

  private final Genomics genomics;
  private final String callsetId;

  private final Pager<SearchVariantsRequest, SearchVariantsResponse> pager =
      new Pager<SearchVariantsRequest, SearchVariantsResponse>() {

        @Override protected SearchVariantsRequest
            setPageToken(SearchVariantsRequest request, String pageToken) {
          return request.setPageToken(pageToken);
        }

        @Override protected SearchVariantsResponse send(SearchVariantsRequest request) {
          try {
            return genomics.variants().search(request).execute();
          } catch (IOException e) {
            throw ExceptionWrapper.wrap(e);
          }
        }

        @Override protected String getNextPageToken(SearchVariantsResponse response) {
          return response.getNextPageToken();
        }
      };

  private ApiCallScanner(Genomics genomics, String callsetId) {
    this.genomics = genomics;
    this.callsetId = callsetId;
  }

  @Override
  public <X> X scan(Callback<? extends X> callback) throws IOException {
    try {
      String datasetId = genomics.callsets().get(callsetId).execute().getDatasetId();
      return callback.scan(
          genomics.variants()
              .getSummary()
              .setDatasetId(datasetId)
              .execute()
              .getContigBounds()
              .stream()
              .map(bound -> new SearchVariantsRequest()
                  .setCallsetIds(Collections.singletonList(callsetId))
                  .setContig(bound.getContig())
                  .setDatasetId(datasetId)
                  .setEndPosition(bound.getUpperBound())
                  .setStartPosition(1L))
              .flatMap(pager::fetchAll)
              .flatMap(response -> response.getVariants().stream())
              .<Call>map(variant -> new Call() {

                    private final com.google.api.services.genomics.model.Call
                        call = Iterables.getOnlyElement(variant.getCalls());

                    @Override public List<String> alternates() {
                      return variant.getAlternateBases();
                    }

                    @Override public String contig() {
                      return variant.getContig();
                    }

                    @Override public boolean equals(Object obj) {
                      return HASH_CODE_AND_EQUALS.equals(this, obj);
                    }

                    @Override public List<Integer> genotype() {
                      return call.getGenotype();
                    }

                    @Override public int hashCode() {
                      return HASH_CODE_AND_EQUALS.hashCode(this);
                    }

                    @Override public Optional<Phaseset> phaseset() {
                      return Optional.ofNullable(call.getPhaseset()).map(ps -> {
                            try {
                              return Call.Phaseset.create(Integer.parseInt(ps));
                            } catch (NumberFormatException e) {
                              return Call.Phaseset.DEFAULT;
                            }
                          });
                    }

                    @Override public int position() {
                      return variant.getPosition().intValue();
                    }

                    @Override public String reference() {
                      return variant.getReferenceBases();
                    }
                  }));
    } catch (ExceptionWrapper e) {
      throw (IOException) e.getCause();
    }
  }
}

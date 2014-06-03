package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.googleapis.batch.BatchRequest;
import com.google.api.client.googleapis.batch.json.JsonBatchCallback;
import com.google.api.client.googleapis.json.GoogleJsonError;
import com.google.api.client.http.HttpHeaders;
import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.collect.AbstractSequentialIterator;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class VariantFetcher {

  public interface Callback<X> {
    X accept(Map<String, Callset> callsets, FluentIterable<Variant> variants) throws Exception;
  }

  private enum VariantFetcherStrategy {

    THE_RIGHT_WAY {

      @Override
      Map<String, Callset> getCallsets(Genomics genomics, Iterable<String> callsetIds)
          throws IOException {
        BatchRequest batch = genomics.batch();
        final ImmutableMap.Builder<String, Callset> callsets = ImmutableMap.builder();
        for (final String callsetId : callsetIds) {
          genomics.callsets().get(callsetId).queue(batch,
              new JsonBatchCallback<Callset>() {

                @Override
                public void onSuccess(Callset callset, HttpHeaders httpHeaders) {
                  callsets.put(callsetId, callset);
                }

                @Override
                public void onFailure(GoogleJsonError error, HttpHeaders httpHeaders)
                    throws IOException {
                  throw new IOException(error.toPrettyString());
                }
              });
        }
        batch.execute();
        return callsets.build();
      }

      @Override
      Map<String, Long> getContigs(Genomics genomics, String datasetId) throws IOException {
        ImmutableMap.Builder<String, Long> contigs = ImmutableMap.builder();
        for (ContigBound contigBound : genomics.variants()
            .getSummary()
            .setDatasetId(datasetId)
            .execute()
            .getContigBounds()) {
          contigs.put(contigBound.getContig(), contigBound.getUpperBound());
        }
        return contigs.build();
      }
    },

    TEMPORARY_HACK {

      @Override
      Map<String, Callset> getCallsets(Genomics genomics, Iterable<String> callsetIds)
          throws IOException {
        ImmutableMap.Builder<String, Callset> callsets = ImmutableMap.builder();
        Set<String> callsetIdSet = ImmutableSet.copyOf(callsetIds);
        for (Callset callset : genomics.callsets()
            .search(new SearchCallsetsRequest()
                .setDatasetIds(Collections.singletonList("208145944813")))
            .execute()
            .getCallsets()) {
          String callsetId = callset.getId();
          if (callsetIdSet.contains(callsetId)) {
            callsets.put(callsetId, callset);
          }
        }
        return callsets.build();
      }

      @Override
      Map<String, Long> getContigs(Genomics genomics, String datasetId) {
        return ImmutableMap.of("20", 63025520L);
      }
    };

    abstract Map<String, Callset> getCallsets(Genomics genomics, Iterable<String> callsetIds)
        throws IOException;

    abstract Map<String, Long> getContigs(Genomics genomics, String datasetId)
        throws IOException;
  }

  private static class RequestResponsePair {

    static final Function<RequestResponsePair, SearchVariantsResponse> RESPONSE =
        new Function<RequestResponsePair, SearchVariantsResponse>() {
          @Override
          public SearchVariantsResponse apply(RequestResponsePair pair) {
            return pair.response;
          }
        };

    static RequestResponsePair of(SearchVariantsRequest request) {
      return of(request, null);
    }

    static RequestResponsePair of(SearchVariantsRequest request, SearchVariantsResponse response) {
      return new RequestResponsePair(request, response);
    }

    private final SearchVariantsRequest request;
    private final SearchVariantsResponse response;

    private RequestResponsePair(SearchVariantsRequest request, SearchVariantsResponse response) {
      this.request = request;
      this.response = response;
    }

    Optional<SearchVariantsRequest> request() {
      return Optional.fromNullable(request);
    }
  }

  private static class RuntimeIOException extends RuntimeException {

    static RuntimeIOException create(IOException cause) {
      return new RuntimeIOException(cause);
    }

    private RuntimeIOException(IOException cause) {
      super(cause);
    }

    @Override
    public synchronized IOException getCause() {
      return (IOException) super.getCause();
    }
  }

  private static final Function<SearchVariantsResponse, Iterable<Variant>>
      SEARCH_VARIANTS_RESPONSE_GET_VARIANTS =
          new Function<SearchVariantsResponse, Iterable<Variant>>() {
            @Override public Iterable<Variant> apply(SearchVariantsResponse response) {
              return response.getVariants();
            }
          };

  public static VariantFetcher create(Genomics genomics) {
    return new VariantFetcher(genomics, genomics.variants());
  }

  private final Genomics genomics;
  private final Genomics.Variants variants;

  private VariantFetcher(Genomics genomics, Genomics.Variants variants) {
    this.genomics = genomics;
    this.variants = variants;
  }

  public <X> X fetchVariants(final List<String> callsetIds, Callback<? extends X> callback)
      throws Exception {
    try {
      VariantFetcherStrategy strategy = VariantFetcherStrategy.TEMPORARY_HACK;
      Map<String, Callset> callsets = strategy.getCallsets(genomics, callsetIds);
      final String datasetId = getDatasetId(callsets);
      return callback.accept(callsets, FluentIterable
          .from(strategy.getContigs(genomics, datasetId).entrySet())
          .transformAndConcat(
              new Function<Map.Entry<String, Long>, Iterable<RequestResponsePair>>() {
                @Override public Iterable<RequestResponsePair> apply(
                    final Map.Entry<String, Long> entry) {
                  return new Iterable<RequestResponsePair>() {
                    @Override public Iterator<RequestResponsePair> iterator() {
                      return new AbstractSequentialIterator<RequestResponsePair>(
                          RequestResponsePair.of(new SearchVariantsRequest()
                              .setDatasetId(datasetId)
                              .setContig(entry.getKey())
                              .setCallsetIds(callsetIds)
                              .setStartPosition(1L)
                              .setEndPosition(entry.getValue()))) {
                        @Override
                        protected RequestResponsePair computeNext(RequestResponsePair previous) {
                          return previous.request()
                              .transform(
                                  new Function<SearchVariantsRequest, RequestResponsePair>() {
                                    @Override public RequestResponsePair apply(
                                        final SearchVariantsRequest request) {
                                      try {
                                        SearchVariantsResponse response =
                                            variants.search(request).execute();
                                        return RequestResponsePair.of(
                                            Optional.fromNullable(response.getNextPageToken())
                                                .transform(
                                                    new Function<String, SearchVariantsRequest>() {
                                                      @Override public SearchVariantsRequest apply(
                                                          String pageToken) {
                                                        return request.setPageToken(pageToken);
                                                      }
                                                    })
                                                .orNull(),
                                            response);
                                      } catch (IOException e) {
                                        throw RuntimeIOException.create(e);
                                      }
                                    }
                                  })
                              .orNull();
                        }
                      };
                    }
                  };
                }
              })
          .transform(RequestResponsePair.RESPONSE)
          .filter(Predicates.notNull())
          .transformAndConcat(SEARCH_VARIANTS_RESPONSE_GET_VARIANTS));
    } catch (RuntimeIOException e) {
      throw e.getCause();
    }
  }

  private static String getDatasetId(Map<?, Callset> callsets) {
    ImmutableMultimap.Builder<String, String> builder = ImmutableMultimap.builder();
    for (Callset callset : callsets.values()) {
      builder.put(callset.getDatasetId(), callset.getId());
    }
    Map<String, Collection<String>> callsetIdsByDatasetId = builder.build().asMap();
    switch (callsetIdsByDatasetId.size()) {
      case 0:
        throw new IllegalArgumentException(
            "At least one callset must be given in order to determine the datasetId");
      case 1:
        return Iterables.getOnlyElement(callsetIdsByDatasetId.keySet());
      default:
        throw new IllegalArgumentException(String.format(
            "All callsets must belong to the same dataset, but %s",
            Joiner.on(", ").join(Iterables.transform(
                callsetIdsByDatasetId.entrySet(),
                new Function<Map.Entry<String, Collection<String>>, String>() {
                  @Override public String apply(Map.Entry<String, Collection<String>> entry) {
                    return String.format(
                        "callsets { %s } belong to dataset %s",
                        Joiner.on(", ").join(entry.getValue()),
                        entry.getKey());
                  }
                }))));
    }
  }
}

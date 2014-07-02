package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.googleapis.batch.BatchRequest;
import com.google.api.client.googleapis.batch.json.JsonBatchCallback;
import com.google.api.client.googleapis.json.GoogleJsonError;
import com.google.api.client.http.HttpHeaders;
import com.google.api.client.repackaged.com.google.common.base.Joiner;
import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;

import com.beust.jcommander.internal.Lists;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class VariantFetcher {

  public interface Callback<X> {
    X accept(Map<String, Callset> callsets, FluentIterable<Variant> variants) throws Exception;
  }

  private static class ExceptionWrapper extends RuntimeException {

    ExceptionWrapper(IOException cause) {
      super(cause);
    }

    @Override public IOException getCause() {
      return (IOException) super.getCause();
    }
  }

  private static final Comparator<ContigBound> CONTIG_BOUND_COMPARATOR = Ordering.natural()
      .onResultOf(
          new Function<ContigBound, Long>() {
            @Override public Long apply(ContigBound bound) {
              return bound.getUpperBound();
            }
          })
      .reverse();

  private static final Function<SearchVariantsResponse, Iterable<Variant>> GET_VARIANTS =
      new Function<SearchVariantsResponse, Iterable<Variant>>() {
        @Override public Iterable<Variant> apply(SearchVariantsResponse response) {
          return response.getVariants();
        }
      };

  private static final Function<Iterable<Callset>, Map<String, Callset>> INDEX_BY_ID =
      new Function<Iterable<Callset>, Map<String, Callset>>() {

        private final Function<Callset, String> getId =
            new Function<Callset, String>() {
              @Override public String apply(Callset callset) {
                return callset.getId();
              }
            };

        @Override public Map<String, Callset> apply(Iterable<Callset> callsets) {
          return FluentIterable.from(callsets).uniqueIndex(getId);
        }
      };

  public static VariantFetcher create(Genomics genomics) {
    return create(genomics, Optional.<List<Callset>>absent());
  }

  public static VariantFetcher create(Genomics genomics, Callset... callsets) {
    return create(genomics, Optional.of(Arrays.asList(callsets)));
  }

  private static VariantFetcher create(Genomics genomics, Optional<List<Callset>> callsets) {
    final Genomics.Variants variants = genomics.variants();
    return new VariantFetcher(
        genomics,
        variants,
        new Pager<SearchVariantsRequest, SearchVariantsResponse>() {

          @Override protected String getNextPageToken(SearchVariantsResponse response) {
            return response.getNextPageToken();
          }

          @Override protected SearchVariantsResponse send(SearchVariantsRequest request) {
            try {
              return variants.search(request).execute();
            } catch (IOException e) {
              throw new ExceptionWrapper(e);
            }
          }

          @Override protected SearchVariantsRequest setPageToken(
              SearchVariantsRequest request, String pageToken) {
            return request.setPageToken(pageToken);
          }
        },
        callsets.transform(INDEX_BY_ID));
  }

  private static String getDatasetId(Map<String, Callset> callsets) {
    ImmutableMultimap.Builder<String, String> builder = ImmutableMultimap.builder();
    for (Map.Entry<String, Callset> entry : callsets.entrySet()) {
      builder.put(entry.getValue().getDatasetId(), entry.getKey());
    }
    Map<String, Collection<String>> callsetIdsByDatasetId = builder.build().asMap();
    switch (callsetIdsByDatasetId.size()) {
      case 0:
        throw new IllegalArgumentException("No callsets to extract datasetId from");
      case 1:
        return Iterables.getOnlyElement(callsetIdsByDatasetId.keySet());
      default:
        throw new IllegalArgumentException(String.format(
            "Not all callsets belong to the same dataset: %s.",
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
  private final Optional<Map<String, Callset>> callsets;
  private final Genomics genomics;
  private final Pager<SearchVariantsRequest, SearchVariantsResponse> searchVariants;

  private final Genomics.Variants variants;

  private VariantFetcher(
      Genomics genomics,
      Genomics.Variants variants,
      Pager<SearchVariantsRequest, SearchVariantsResponse> searchVariants,
      Optional<Map<String, Callset>> callsets) {
    this.genomics = genomics;
    this.variants = variants;
    this.searchVariants = searchVariants;
    this.callsets = callsets;
  }

  public <X> X fetchVariants(final List<String> callsetIds, Callback<? extends X> callback)
      throws Exception {
    try {
      Map<String, Callset> callsets = this.callsets.isPresent()
          ? Maps.filterKeys(this.callsets.get(), Predicates.in(callsetIds))
          : getCallsets(callsetIds);
      final String datasetId = getDatasetId(callsets);
      return callback.accept(
          callsets,
          FluentIterable.from(getContigs(datasetId).entrySet())
              .transform(
                  new Function<Map.Entry<String, Long>, SearchVariantsRequest>() {
                    @Override public SearchVariantsRequest apply(Map.Entry<String, Long> entry) {
                      return new SearchVariantsRequest()
                          .setDatasetId(datasetId)
                          .setCallsetIds(callsetIds)
                          .setContig(entry.getKey())
                          .setStartPosition(1L)
                          .setEndPosition(entry.getValue());
                    }
                  })
              .transformAndConcat(searchVariants)
              .transformAndConcat(GET_VARIANTS));
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }

  private Map<String, Callset> getCallsets(List<String> callsetIds) throws IOException {
    BatchRequest batch = genomics.batch();
    final ImmutableMap.Builder<String, Callset> callsets = ImmutableMap.builder();
    for (final String callsetId : callsetIds) {
      genomics.callsets().get(callsetId).queue(batch,
          new JsonBatchCallback<Callset>() {

            @Override public void onFailure(GoogleJsonError error, HttpHeaders headers)
                throws IOException {
              throw new IOException(error.toPrettyString());
            }

            @Override public void onSuccess(Callset callset, HttpHeaders headers) {
              callsets.put(callsetId, callset);
            }
          });
    }
    batch.execute();
    return callsets.build();
  }

  private Map<String, Long> getContigs(String datasetId) throws IOException {
    List<ContigBound> bounds = Lists.newArrayList(
        variants.getSummary().setDatasetId(datasetId).execute().getContigBounds());
    Collections.sort(bounds, CONTIG_BOUND_COMPARATOR);
    ImmutableMap.Builder<String, Long> contigs = ImmutableMap.builder();
    for (ContigBound bound : bounds) {
      contigs.put(bound.getContig(), bound.getUpperBound());
    }
    return contigs.build();
  }
}

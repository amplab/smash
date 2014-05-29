package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.googleapis.batch.BatchRequest;
import com.google.api.client.googleapis.batch.json.JsonBatchCallback;
import com.google.api.client.googleapis.json.GoogleJsonError;
import com.google.api.client.http.HttpHeaders;
import com.google.api.client.repackaged.com.google.common.base.Joiner;
import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.Iterables;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class VariantFetcher {

  public interface Callback<X> {
    X accumulate(X accumulator, Variant variant) throws IOException;
  }

  public static VariantFetcher create(Genomics genomics) {
    return new VariantFetcher(genomics, genomics.callsets(), genomics.variants());
  }

  private final Genomics genomics;
  private final Genomics.Callsets callsets;
  private final Genomics.Variants variants;

  private VariantFetcher(
      Genomics genomics,
      Genomics.Callsets callsets,
      Genomics.Variants variants) {
    this.genomics = genomics;
    this.callsets = callsets;
    this.variants = variants;
  }

  public <X> X fetchVariants(
      SearchVariantsRequest prototype,
      List<String> contigs,
      X initial,
      Callback<X> callback) throws IOException {
    X accumulator = initial;
    for (String contig : contigs) {
      SearchVariantsRequest request = clone(prototype).setContig(contig);
      for (boolean condition = true; condition; ) {
        SearchVariantsResponse response = variants.search(request).execute();
        for (Variant variant : Optional.fromNullable(response.getVariants())
            .or(Collections.<Variant>emptyList())) {
          accumulator = callback.accumulate(accumulator, variant);
        }
        Optional<String> pageToken = Optional.fromNullable(response.getNextPageToken());
        if (condition = pageToken.isPresent()) {
          request.setPageToken(pageToken.get());
        }
      }
    }
    return accumulator;
  }

  public String getDatasetId(List<String> callsetIds) throws IOException {
    Preconditions.checkArgument(!callsetIds.isEmpty(), "callsetIds was empty");
    BatchRequest batch = genomics.batch();
    final ImmutableMultimap.Builder<String, String> builder = ImmutableMultimap.builder();
    for (String callsetId : callsetIds) {
      callsets.get(callsetId).queue(batch,
          new JsonBatchCallback<Callset>() {

            @Override public void onFailure(GoogleJsonError error, HttpHeaders httpHeaders)
                throws IOException {
              throw new IOException(error.toPrettyString());
            }

            @Override public void onSuccess(Callset callset, HttpHeaders httpHeaders) {
              builder.put(callset.getDatasetId(), callset.getId());
            }
          });
    }
    batch.execute();
    Map<String, Collection<String>> callsetsByDatasetId = builder.build().asMap();
    if (1 < callsetIds.size()) {
      throw new IllegalStateException(String.format(
          "Every callset specified must belong to the same dataset, but %s",
          Joiner.on(", ").join(Iterables.transform(
              callsetsByDatasetId.entrySet(),
              new Function<Map.Entry<String, Collection<String>>, String>() {
                @Override public String apply(Map.Entry<String, Collection<String>> entry) {
                  return String.format(
                      "callsets { %s } belong to dataset %s",
                      Joiner.on(", ").join(entry.getValue()),
                      entry.getKey());
                }
              }))));
    }
    return Iterables.getOnlyElement(callsetsByDatasetId.keySet());
  }

  /**
   * Because {@code SearchVariantsRequest.clone()} is busted.
   */
  static SearchVariantsRequest clone(SearchVariantsRequest prototype) {
    SearchVariantsRequest copy = new SearchVariantsRequest();
    Optional<List<String>> callsetIds = Optional.fromNullable(prototype.getCallsetIds());
    if (callsetIds.isPresent()) {
      copy.setCallsetIds(callsetIds.get());
    }
    Optional<List<String>> callsetNames = Optional.fromNullable(prototype.getCallsetNames());
    if (callsetNames.isPresent()) {
      copy.setCallsetNames(callsetNames.get());
    }
    Optional<String> contig = Optional.fromNullable(prototype.getContig());
    if (contig.isPresent()) {
      copy.setContig(contig.get());
    }
    Optional<String> datasetId = Optional.fromNullable(prototype.getDatasetId());
    if (datasetId.isPresent()) {
      copy.setDatasetId(datasetId.get());
    }
    Optional<Long> endPosition = Optional.fromNullable(prototype.getEndPosition());
    if (endPosition.isPresent()) {
      copy.setEndPosition(endPosition.get());
    }
    Optional<BigInteger> maxResults = Optional.fromNullable(prototype.getMaxResults());
    if (maxResults.isPresent()) {
      copy.setMaxResults(maxResults.get());
    }
    Optional<String> pageToken = Optional.fromNullable(prototype.getPageToken());
    if (pageToken.isPresent()) {
      copy.setPageToken(pageToken.get());
    }
    Optional<Long> startPosition = Optional.fromNullable(prototype.getStartPosition());
    if (startPosition.isPresent()) {
      copy.setStartPosition(startPosition.get());
    }
    Optional<String> variantName = Optional.fromNullable(prototype.getVariantName());
    if (variantName.isPresent()) {
      copy.setVariantName(variantName.get());
    }
    return copy;
  }
}

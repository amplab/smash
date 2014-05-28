package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Optional;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Collections;
import java.util.List;

public class VariantFetcher {

  public interface Callback<X> {
    X accumulate(X accumulator, Variant variant) throws IOException;
  }

  public static VariantFetcher create(Genomics genomics) {
    return new VariantFetcher(genomics.variants());
  }

  private final Genomics.Variants variants;

  private VariantFetcher(Genomics.Variants variants) {
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

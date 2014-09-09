/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.berkeley.cs.amplab.calldiff;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.cloud.genomics.utils.Paginator;
import com.google.common.collect.Iterables;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.StreamSupport;

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
  private final Paginator.Variants searchVariants;

  private ApiCallScanner(Genomics genomics, String callsetId) {
    this.genomics = genomics;
    this.callsetId = callsetId;
    this.searchVariants = Paginator.Variants.create(genomics);
  }

  @Override
  public <X> X scan(Callback<? extends X> callback) throws IOException {
    try {
      String variantsetId = genomics.callsets().get(callsetId).execute().getVariantsetId();
      return callback.scan(
          genomics.variants()
              .getSummary()
              .setVariantsetId(variantsetId)
              .execute()
              .getContigBounds()
              .stream()
              .map(bound -> new SearchVariantsRequest()
                  .setCallsetIds(Collections.singletonList(callsetId))
                  .setContig(bound.getContig())
                  .setVariantsetId(variantsetId)
                  .setEndPosition(bound.getUpperBound())
                  .setStartPosition(1L))
              .flatMap(request -> StreamSupport.stream(
                  Spliterators.spliteratorUnknownSize(
                      searchVariants.search(request, search -> {}).iterator(),
                      Spliterator.DISTINCT | Spliterator.IMMUTABLE | Spliterator.NONNULL),
                  false))
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

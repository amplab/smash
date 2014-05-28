package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.runners.MockitoJUnitRunner;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

@RunWith(MockitoJUnitRunner.class)
public class VariantFetcherTest {

  private static final Function<Variant, String> GET_CONTIG =
      new Function<Variant, String>() {
        @Override public String apply(Variant variant) {
          return variant.getContig();
        }
      };

  @Mock private Genomics.Variants variants;
  private VariantFetcher fetcher;

  @Before
  public void setUp() {
    Genomics genomics = Mockito.mock(Genomics.class);
    when(genomics.variants()).thenReturn(variants);
    fetcher = VariantFetcher.create(genomics);
  }

  @Test
  public void testFetchVariants() throws IOException {
    int maxResults = 2;
    SearchVariantsRequest request = new SearchVariantsRequest()
        .setMaxResults(new BigInteger(String.valueOf(maxResults)));
    List<String> contigs = Arrays.asList("1", "2");
    List<Variant> expected = Arrays.asList(
        variant(1, 1),
        variant(2, 1),
        variant(3, 1),
        variant(4, 1),
        variant(5, 1),
        variant(6, 2),
        variant(7, 2),
        variant(8, 2));
    for (String contig : contigs) {
      Optional<String> pageToken = Optional.absent();
      for (
          Iterator<List<Variant>> iterator = Iterables
              .partition(
                  Iterables.filter(
                      expected,
                      Predicates.compose(Predicates.equalTo(contig), GET_CONTIG)),
                  maxResults)
              .iterator();
          iterator.hasNext();) {
        SearchVariantsRequest request1 = VariantFetcher.clone(request)
            .setContig(contig);
        if (pageToken.isPresent()) {
          request1.setPageToken(pageToken.get());
        }
        SearchVariantsResponse response = new SearchVariantsResponse()
            .setVariants(iterator.next());
        if (iterator.hasNext()) {
          response.setNextPageToken((pageToken = Optional.of(pageToken.isPresent()
              ? String.valueOf(Integer.parseInt(pageToken.get()) + 1)
              : "2")).get());
        }
        Genomics.Variants.Search search = mock(Genomics.Variants.Search.class);
        when(search.execute()).thenReturn(response);
        when(this.variants.search(request1)).thenReturn(search);
      }
    }
    assertEquals(expected, fetcher
        .fetchVariants(request, contigs, ImmutableList.<Variant>builder(),
            new VariantFetcher.Callback<ImmutableList.Builder<Variant>>() {
              @Override public ImmutableList.Builder<Variant> accumulate(
                  ImmutableList.Builder<Variant> variants, Variant variant) {
                return variants.add(variant);
              }
            })
        .build());
  }

  private static Variant variant(int id, int contig) {
    return new Variant().setId(String.valueOf(id)).setContig(String.valueOf(contig));
  }
}

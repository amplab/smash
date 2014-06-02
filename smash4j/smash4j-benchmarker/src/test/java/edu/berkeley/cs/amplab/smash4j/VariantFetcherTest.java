package edu.berkeley.cs.amplab.smash4j;

import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchCallsetsResponse;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;
import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Iterables;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.Mock;
import org.mockito.runners.MockitoJUnitRunner;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

@RunWith(MockitoJUnitRunner.class)
public class VariantFetcherTest {

  @Mock private Genomics genomics;
  @Mock private Genomics.Variants variants;

  private VariantFetcher fetcher;

  @Before
  public void setUp() {
    when(genomics.variants()).thenReturn(variants);
    fetcher = VariantFetcher.create(genomics);
  }

  @Test
  public void testFetchVariants() throws Exception {
    class Actual {

      final Map<String, Callset> callsets;
      final List<Variant> variants;

      Actual(Map<String, Callset> callsets, List<Variant> variants) {
        this.callsets = callsets;
        this.variants = variants;
      }
    }
    Map<String, Callset> expectedCallsets = ImmutableMap.of(
        "208145944813-0", new Callset().setId("208145944813-0").setDatasetId("208145944813"),
        "208145944813-1", new Callset().setId("208145944813-1").setDatasetId("208145944813"),
        "208145944813-2", new Callset().setId("208145944813-2").setDatasetId("208145944813"));
    List<Variant> expectedVariants = Arrays.asList(
        variant(1),
        variant(2),
        variant(3),
        variant(4),
        variant(5),
        variant(6),
        variant(7),
        variant(8));
    Genomics.Callsets callsets = mock(Genomics.Callsets.class);
    Genomics.Callsets.Search callsetsSearch = mock(Genomics.Callsets.Search.class);
    when(genomics.callsets()).thenReturn(callsets);
    when(callsets.search(new SearchCallsetsRequest().setDatasetIds(
        Collections.singletonList("208145944813")))).thenReturn(callsetsSearch);
    when(callsetsSearch.execute()).thenReturn(new SearchCallsetsResponse()
        .setCallsets(ImmutableList.copyOf(expectedCallsets.values())));
    Optional<String> pageToken = Optional.absent();
    for (Iterator<List<Variant>> iterator = Iterables.partition(expectedVariants, 3).iterator();
         iterator.hasNext();) {
      List<Variant> page = iterator.next();
      SearchVariantsRequest request = new SearchVariantsRequest()
          .setDatasetId("208145944813")
          .setContig("20")
          .setCallsetIds(Arrays.asList("208145944813-0", "208145944813-1", "208145944813-2"))
          .setStartPosition(1L)
          .setEndPosition(63025520L);
      if (pageToken.isPresent()) {
        request.setPageToken(pageToken.get());
      }
      SearchVariantsResponse response = new SearchVariantsResponse()
          .setVariants(page);
      if (iterator.hasNext()) {
        response.setNextPageToken((pageToken = Optional.of(String.valueOf(
            pageToken.isPresent() ? Integer.parseInt(pageToken.get() + 1) : 2))).get());
      }
      Genomics.Variants.Search variantsSearch = mock(Genomics.Variants.Search.class);
      when(variants.search(request)).thenReturn(variantsSearch);
      when(variantsSearch.execute()).thenReturn(response);
    }
    Actual actual = fetcher.fetchVariants(
        Arrays.asList("208145944813-0", "208145944813-1", "208145944813-2"),
        new VariantFetcher.Callback<Actual>() {
          @Override
          public Actual accept(
              Map<String, Callset> callsets, FluentIterable<Variant> variants) {
            return new Actual(callsets, variants.toList());
          }
        });
    assertEquals(expectedCallsets, actual.callsets);
    assertEquals(expectedVariants, actual.variants);
  }

  private static Variant variant(int id) {
    return new Variant().setId(String.valueOf(id));
  }
}

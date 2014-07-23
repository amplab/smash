package edu.berkeley.cs.amplab.calldiff;

import com.google.common.collect.AbstractSequentialIterator;

import java.util.Optional;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Generic code that understands the {@code pageToken} - {@code nextPageToken} pagination protocol
 * that the genomics API uses.
 */
public abstract class Pager<Q, A> {

  private class Pair {

    final Q request;
    final A response;

    Pair(Q request, A response) {
      this.request = request;
      this.response = response;
    }
  }

  public final Stream<A> fetchAll(Q initial) {
    return StreamSupport
        .stream(
            Spliterators.spliteratorUnknownSize(
                new AbstractSequentialIterator<Pair>(new Pair(initial, null)) {
                  @Override protected Pair computeNext(Pair previous) {
                    return Optional.ofNullable(previous.request)
                        .map(request -> {
                              A response = send(request);
                              return new Pair(
                                  Optional.ofNullable(getNextPageToken(response))
                                      .map(pageToken -> setPageToken(request, pageToken))
                                      .orElse(null),
                                  response);
                            })
                        .orElse(null);
                  }
                },
                Spliterator.DISTINCT | Spliterator.IMMUTABLE | Spliterator.NONNULL),
            false)
        .map(pair -> pair.response)
        .skip(1);
  }

  protected abstract String getNextPageToken(A response);

  protected abstract A send(Q request);

  protected abstract Q setPageToken(Q request, String pageToken);
}

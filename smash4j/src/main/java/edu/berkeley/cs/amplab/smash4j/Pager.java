package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.collect.AbstractSequentialIterator;
import com.google.common.collect.FluentIterable;

import java.util.Iterator;

public abstract class Pager<Request, Response>
    implements Function<Request, FluentIterable<Response>> {

  private class Pair {

    final Request request;
    final Response response;

    Pair(Request request, Response response) {
      this.request = request;
      this.response = response;
    }
  }

  private final Function<Pair, Response> getResponse =
      new Function<Pair, Response>() {
        @Override public Response apply(Pair pair) {
          return pair.response;
        }
      };

  @Override
  public final FluentIterable<Response> apply(final Request initial) {
    return new FluentIterable<Pair>() {
          @Override public Iterator<Pair> iterator() {
            return new AbstractSequentialIterator<Pair>(new Pair(initial, null)) {
                  @Override protected Pair computeNext(Pair previous) {
                    return Optional.fromNullable(previous.request)
                        .transform(
                            new Function<Request, Pair>() {
                              @Override public Pair apply(final Request request) {
                                Response response = send(request);
                                return new Pair(
                                    Optional.fromNullable(getNextPageToken(response))
                                        .transform(new Function<String, Request>() {
                                              @Override public Request apply(String pageToken) {
                                                return setPageToken(request, pageToken);
                                              }
                                            })
                                        .orNull(),
                                    response);
                              }
                            })
                        .orNull();
                  }
            };
          }
        }.transform(getResponse).skip(1);
  }

  protected abstract Response send(Request request);

  protected abstract String getNextPageToken(Response response);

  protected abstract Request setPageToken(Request request, String pageToken);
}

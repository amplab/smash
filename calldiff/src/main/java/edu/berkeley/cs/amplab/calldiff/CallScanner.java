package edu.berkeley.cs.amplab.calldiff;

import java.io.IOException;
import java.util.stream.Stream;

/**
 * An interface encapsulating the details of how to consume a set of calls from an input source.
 */
@FunctionalInterface
public interface CallScanner {

  @FunctionalInterface
  interface Callback<X> {

    X scan(Stream<Call> calls);
  }

  <X> X scan(Callback<? extends X> callback) throws IOException;
}

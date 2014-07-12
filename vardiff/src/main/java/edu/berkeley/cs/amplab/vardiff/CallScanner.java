package edu.berkeley.cs.amplab.vardiff;

import java.io.IOException;
import java.util.stream.Stream;

@FunctionalInterface
public interface CallScanner {

  @FunctionalInterface
  interface Callback<X> {

    X scan(Stream<Call> calls);
  }

  <X> X scan(Callback<? extends X> callback) throws IOException;
}

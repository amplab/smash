package edu.berkeley.cs.amplab.calldiff;

import java.util.function.Function;

/**
 * A type used for submarining checked exceptions through methods that don't declare them, for
 * example, {@link Function#apply}. The top level {@link Main#main} method contains a
 * {@code try-catch} block to catch and unwrap these exceptions.
 */
public class ExceptionWrapper extends RuntimeException {

  public static ExceptionWrapper wrap(Exception cause) {
    return new ExceptionWrapper(cause);
  }

  private ExceptionWrapper(Exception cause) {
    super(cause);
  }

  @Override
  public Exception getCause() {
    return (Exception) super.getCause();
  }
}

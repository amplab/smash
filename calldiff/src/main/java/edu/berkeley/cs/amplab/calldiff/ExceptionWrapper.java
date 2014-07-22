package edu.berkeley.cs.amplab.calldiff;

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

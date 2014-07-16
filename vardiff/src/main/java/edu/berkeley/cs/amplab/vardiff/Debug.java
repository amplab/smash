package edu.berkeley.cs.amplab.vardiff;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class Debug {

  private static BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

  public static <X> X debug(String message, X obj) {
    try {
      System.out.format("%s%n%s%n", message, obj);
      in.readLine();
      return obj;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }
}

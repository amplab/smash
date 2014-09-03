package edu.berkeley.cs.amplab.calldiff;

public class CaseNormalizer {

  static String exceptionMessage(int c) {
    return String.format("Illegal nucleobase: %c", c);
  }

  public static char normalizeCase(int c) {
    switch (c) {
      case 'A':
      case 'a':
        return 'A';
      case 'C':
      case 'c':
        return 'C';
      case 'G':
      case 'g':
        return 'G';
      case 'T':
      case 't':
        return 'T';
      default:
        throw new IllegalArgumentException(exceptionMessage(c));
    }
  }

  public static String normalizeCase(String string) {
    return string.codePoints()
        .map(CaseNormalizer::normalizeCase)
        .collect(StringBuilder::new, StringBuilder::appendCodePoint, StringBuilder::append)
        .toString();
  }
}

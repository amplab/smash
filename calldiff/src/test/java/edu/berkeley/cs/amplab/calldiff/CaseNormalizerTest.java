package edu.berkeley.cs.amplab.calldiff;

import static edu.berkeley.cs.amplab.calldiff.CaseNormalizer.exceptionMessage;
import static edu.berkeley.cs.amplab.calldiff.CaseNormalizer.normalizeCase;
import static java.lang.Character.toUpperCase;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CaseNormalizerTest {

  private static final Set<Character> ALLOWED_CHARS = Stream.of('A', 'C', 'G', 'T')
      .flatMap(c -> Stream.of(c, Character.toLowerCase(c)))
      .collect(Collectors.toSet());

  @Test
  public void testNormalizeCase_chars() {
    for (char testChar = 0; testChar < 128; ++testChar) {
      if (ALLOWED_CHARS.contains(testChar)) {
        assertEquals(toUpperCase(testChar), normalizeCase(testChar));
      } else {
        try {
          normalizeCase(testChar);
          fail(Character.toString(testChar));
        } catch (IllegalArgumentException e) {
          assertEquals(exceptionMessage(testChar), e.getMessage());
        }
      }
    }
  }

  @Test
  public void testNormalizeCase_strings() {
    assertEquals("ACGT", normalizeCase("ACGT"));
    assertEquals("ACGT", normalizeCase("acgt"));
    try {
      normalizeCase("foo");
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals(exceptionMessage('f'), e.getMessage());
    }
  }
}

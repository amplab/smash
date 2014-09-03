/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

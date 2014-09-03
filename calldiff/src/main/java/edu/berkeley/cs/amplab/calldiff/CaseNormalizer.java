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

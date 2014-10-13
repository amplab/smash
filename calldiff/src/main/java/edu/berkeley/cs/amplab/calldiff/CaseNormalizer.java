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
      case 'U':
      case 'u':
        return 'U';
      case 'R':
      case 'r':
        return 'R';
      case 'Y':
      case 'y':
        return 'Y';
      case 'K':
      case 'k':
        return 'K';
      case 'M':
      case 'm':
        return 'M';
      case 'S':
      case 's':
        return 'S';
      case 'W':
      case 'w':
        return 'W';
      case 'B':
      case 'b':
        return 'B';
      case 'D':
      case 'd':
        return 'D';
      case 'H':
      case 'h':
        return 'H';
      case 'V':
      case 'v':
        return 'V';
      case 'N':
      case 'n':
        return 'N';
      case 'X':
      case 'x':
        return 'X';
      case '-':
        return '-';
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

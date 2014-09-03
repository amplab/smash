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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;

import org.junit.Test;

import edu.berkeley.cs.amplab.calldiff.FastaReader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

/**
 * Unit test for {@link FastaReader}
 */
public class FastaReaderTest {

  private static class Random {

    static Random create(int lowerBound, int upperBound) {
      return new Random(lowerBound, upperBound);
    }

    private final int lowerBound;
    private final int upperBound;

    private Random(int lowerBound, int upperBound) {
      this.lowerBound = lowerBound;
      this.upperBound = upperBound;
    }

    int nextInt() {
      return RANDOM.nextInt(upperBound - lowerBound + 1) + lowerBound;
    }
  }

  private static final int
      BASES_PER_LINE = 80;
  private static final Random
      CONTIG_LENGTH = Random.create(500, 1000),
      NUMBER_OF_CONTIGS = Random.create(10, 20),
      NUMBER_OF_QUERIES = Random.create(1000, 2000),
      QUERY_LENGTH = Random.create(60, 120);
  private static final java.util.Random
      RANDOM = new java.util.Random();

  /**
   * Create a random FASTA file, perform some random queries against it, and assert that what we
   * get back from the queries was what we expected.
   */
  @Test
  public void testFastaReader() throws Exception {
    final Map<String, String> fasta = randomFasta();
    List<Map.Entry<String, String>> index = ImmutableList.copyOf(fasta.entrySet());
    int indexSize = index.size();
    assertSame(
        this,
        FastaReader.create(saveToTmpFile(fasta))
            .read(fastaFile -> {
              int numberOfQueries = NUMBER_OF_QUERIES.nextInt();
              for (int i = 0; i < numberOfQueries; ++i) {
                Map.Entry<String, String> entry = index.get(RANDOM.nextInt(indexSize));
                String contigName = entry.getKey(), contig = entry.getValue();
                int queryLength = QUERY_LENGTH.nextInt(),
                queryStart = Random.create(0, contig.length() - queryLength - 1).nextInt(),
                queryEnd = queryStart + queryLength;
                assertEquals(
                    fasta.get(contigName).substring(queryStart, queryEnd),
                    fastaFile.get(contigName, queryStart, queryEnd));
              }
              return FastaReaderTest.this;
            }));
  }

  private static Map<String, String> randomFasta() {
    ImmutableMap.Builder<String, String> fasta = ImmutableMap.builder();
    int numberOfContigs = NUMBER_OF_CONTIGS.nextInt();
    for (int i = 0; i < numberOfContigs; ++i) {
      int contigLength = CONTIG_LENGTH.nextInt();
      StringBuilder contig = new StringBuilder();
      for (int j = 0; j < contigLength; ++j) {
        contig.append("ACGT".charAt(RANDOM.nextInt(4)));
      }
      fasta.put(String.format("SEQUENCE_%d", i + 1), contig.toString());
    }
    return fasta.build();
  }

  private static File saveToTmpFile(Map<String, String> fasta) throws IOException {
    File file = File.createTempFile("tmp", ".fa");
    file.deleteOnExit();
    try (PrintStream out = new PrintStream(new FileOutputStream(file))) {
      for (Map.Entry<String, String> entry : fasta.entrySet()) {
        out.format(">%s%n", entry.getKey());
        String contig = entry.getValue();
        int contigLength = contig.length(), fullLines = contigLength / BASES_PER_LINE;
        for (int i = 0; i < fullLines; ++i) {
          out.println(contig.substring(BASES_PER_LINE * i, BASES_PER_LINE * (i + 1)));
        }
        if (0 < contigLength % BASES_PER_LINE) {
          out.println(contig.substring(BASES_PER_LINE * fullLines));
        }
      }
    }
    return file;
  }
}

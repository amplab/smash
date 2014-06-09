package edu.berkeley.cs.amplab.fastaparser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;

import org.junit.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

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

  @Test
  public void testFastaReader() throws Exception {
    final Map<String, String> fasta = randomFasta();
    assertSame(this, FastaReader.create(saveToTmpFile(fasta)).read(
        new FastaReader.Callback<FastaReaderTest>() {

          private final List<Map.Entry<String, String>>
              index = ImmutableList.copyOf(fasta.entrySet());
          private final int indexSize = index.size();

          @Override public FastaReaderTest read(FastaReader.Callback.FastaFile fastaFile)
              throws Exception {
            int numberOfQueries = NUMBER_OF_QUERIES.nextInt();
            for (int i = 0; i < numberOfQueries; ++i) {
              Map.Entry<String, String> entry = index.get(RANDOM.nextInt(indexSize));
              String contigName = entry.getKey(), contig = entry.getValue();
              int queryLength = QUERY_LENGTH.nextInt(),
                  queryStart = Random.create(0, contig.length() - queryLength - 1).nextInt(),
                  queryEnd = queryStart + queryLength;
              Optional<String> optional = fastaFile.get(contigName, queryStart, queryEnd);
              assertTrue(optional.isPresent());
              assertEquals(fasta.get(contigName).substring(queryStart, queryEnd), optional.get());
            }
            return FastaReaderTest.this;
          }
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
        for (int i = 0; i < fullLines;
            out.println(contig.substring(BASES_PER_LINE * i, BASES_PER_LINE * (++i))));
        if (0 < contigLength % BASES_PER_LINE) {
          out.println(contig.substring(BASES_PER_LINE * fullLines));
        }
      }
    }
    return file;
  }
}

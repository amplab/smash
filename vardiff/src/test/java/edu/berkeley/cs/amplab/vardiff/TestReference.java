package edu.berkeley.cs.amplab.vardiff;

import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.google.common.base.Throwables;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;

public class TestReference {

  private static final Supplier<FastaReader> READER = Suppliers.memoize(() -> {
    try {
      File fastaFile = File.createTempFile("temp", ".fasta");
      fastaFile.deleteOnExit();
      try (PrintStream out = new PrintStream(new FileOutputStream(fastaFile))) {
        String dna = "ACGT";
        int dnaLength = dna.length(),
            lineLength = 80;
        Random random = new Random();
        Maps.transformValues(
                ImmutableMap.of("chr1", 100, "chr2", 100),
                length -> {
                  StringBuilder contig = new StringBuilder();
                  for (int i = 0; i < length; ++i) {
                    contig.append(dna.charAt(random.nextInt(dnaLength)));
                  }
                  return contig.toString();
                })
            .entrySet()
            .stream()
            .forEach(entry -> {
              String contig = entry.getValue();
              int contigLength = contig.length(),
                  fullLines = contigLength / lineLength,
                  partialLine = contigLength % lineLength;
              out.format(">%s%n", entry.getKey());
              for (int i = 0; i < fullLines; ++i) {
                out.println(contig.substring(lineLength * i, lineLength * (i + 1)));
              }
              if (0 < partialLine) {
                out.println(contig.substring(lineLength * fullLines, contigLength));
              }
            });
      }
      return FastaReader.create(fastaFile);
    } catch (IOException e) {
      throw Throwables.propagate(e);
    }
  });

  public static FastaReader reader() {
    return READER.get();
  }
}

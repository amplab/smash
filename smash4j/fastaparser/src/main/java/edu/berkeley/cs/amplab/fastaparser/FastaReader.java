package edu.berkeley.cs.amplab.fastaparser;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Map;

public class FastaReader {

  public interface Callback<X> {

    public interface FastaFile {

      public enum Orientation {

        FORWARD {
          @Override public String apply(String sequence) {
            return sequence;
          }
        },

        REVERSE {

          @Override public String apply(String sequence) {
            StringBuilder builder = new StringBuilder();
            for (
                int i = sequence.length() - 1;
                0 <= i;
                builder.append(compliment(sequence.charAt(i--))));
            return builder.toString();
          }

          private char compliment(char c) {
            switch (c) {
              case 'A': return 'T';
              case 'a': return 't';
              case 'C': return 'G';
              case 'c': return 'g';
              case 'T': return 'A';
              case 't': return 'a';
              case 'G': return 'C';
              case 'g': return 'c';
              default : return c;
            }
          }
        };

        abstract String apply(String sequence);
      }

      String get(String contigName, int beginIndex, int endIndex, Orientation orientation);
    }

    X read(Map<String, Integer> info, FastaFile fastaFile) throws Exception;
  }

  public static FastaReader create(File fastaFile) throws IOException {
    File faiFile = fastaFile.toPath()
        .getParent()
        .resolve(String.format("%s.fai", fastaFile.getName()))
        .toFile();
    return faiFile.exists() && faiFile.isFile() && faiFile.canRead()
        ? create(fastaFile, faiFile)
        : new FastaReader(fastaFile, FastaIndex.create(fastaFile));
  }

  public static FastaReader create(File fastaFile, File faiFile) throws IOException {
    return new FastaReader(fastaFile, FastaIndex.read(faiFile));
  }

  private final File fastaFile;
  private final FastaIndex index;

  private FastaReader(File fastaFile, FastaIndex index) {
    this.fastaFile = fastaFile;
    this.index = index;
  }

  public <X> X read(Callback<? extends X> callback) throws Exception {
    try (RandomAccessFile file = new RandomAccessFile(fastaFile, "r")) {
      try (FileChannel channel = file.getChannel()) {
        class Contig {

          private final ByteBuffer contig;
          private final String contigName;
          private final int basesPerLine;

          Contig(ByteBuffer contig, String contigName, int basesPerLine) {
            this.contig = contig;
            this.contigName = contigName;
            this.basesPerLine = basesPerLine;
          }

          String get(int beginIndex, int endIndex) {
            StringBuilder builder = new StringBuilder();
            for (int end = accountForNewlines(endIndex), i = accountForNewlines(beginIndex);
                i < end; ++i) {
              char c = (char) contig.get(i);
              if (Character.isAlphabetic(c) || '-' == c) {
                builder.append(c);
              } else if ('\n' != c) {
                throw new IllegalStateException(String.format(
                    "Illegal character %c on contig \"%s\" at position %d", c, contigName, i));
              }
            }
            return builder.toString();
          }

          int length() {
            int limit = contig.limit();
            return limit - limit / basesPerLine;
          }

          private int accountForNewlines(int i) {
            return i + i / basesPerLine;
          }
        }
        ImmutableMap.Builder<String, Contig> builder = ImmutableMap.builder();
        for (FastaIndex.Entry entry : index.entries()) {
          long length = entry.length();
          int bases = entry.bases();
          String name = entry.name();
          builder.put(name, new Contig(channel.map(FileChannel.MapMode.READ_ONLY, entry.offset(),
              (length / bases) * entry.bytes() + length % bases), name, entry.bases()));
        }
        final Map<String, Contig> chromosomes = builder.build();
        return callback.read(
            Maps.transformValues(chromosomes,
                new Function<Contig, Integer>() {
                  @Override public Integer apply(Contig contig) {
                    return contig.length();
                  }
                }),
            new Callback.FastaFile() {
              @Override public String get(String contigName, int beginIndex, int endIndex,
                  Callback.FastaFile.Orientation orientation) {
                return orientation.apply(chromosomes.get(contigName).get(beginIndex, endIndex));
              }
            });
      }
    }
  }
}

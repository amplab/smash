package edu.berkeley.cs.amplab.smash4j.fasta;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Map;

public class FastaReader {

  public interface Callback<X> {
    X read(Map<String, Integer> info, FastaFile fastaFile) throws Exception;
  }

  public interface FastaFile {
    String get(String contigName, int beginIndex, int endIndex);
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

          private final int basesPerLine;
          private final ByteBuffer contig;
          private final String contigName;
          private final int maxIndex;

          Contig(ByteBuffer contig, int contigSize, String contigName, int basesPerLine) {
            this.contig = contig;
            this.maxIndex = contigSize;
            this.contigName = contigName;
            this.basesPerLine = basesPerLine;
          }

          private int accountForNewlines(int i) {
            return i + i / basesPerLine;
          }

          String get(int beginIndex, int endIndex) {
            StringBuilder builder = new StringBuilder();
            for (int end = Math.min(maxIndex, accountForNewlines(endIndex)),
                i = Math.max(0, accountForNewlines(beginIndex)); i < end; ++i) {
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
        }
        ImmutableMap.Builder<String, Contig> builder = ImmutableMap.builder();
        for (FastaIndex.Entry entry : index.entries()) {
          long length = entry.length();
          int bases = entry.bases();
          String name = entry.name();
          long position = entry.offset();
          MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_ONLY, position, Math
              .min((length / bases) * entry.bytes() + length % bases, file.length() - position));
          builder.put(name, new Contig(buffer, buffer.limit() - 1, name, entry.bases()));
        }
        final Map<String, Contig> chromosomes = builder.build();
        return callback.read(
            Maps.transformValues(chromosomes,
                new Function<Contig, Integer>() {
                  @Override public Integer apply(Contig contig) {
                    return contig.length();
                  }
                }),
            new FastaFile() {
              @Override public String get(String contigName, int beginIndex, int endIndex) {
                return chromosomes.get(contigName).get(beginIndex, endIndex);
              }
            });
      }
    }
  }
}

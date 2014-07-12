package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.ImmutableMap;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class FastaReader {

  @FunctionalInterface
  public interface Callback<X> {
    X read(FastaFile fastaFile);
  }

  public interface FastaFile {

    int contigLength(String contig);

    Set<String> contigs();

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

  public <X> X read(Callback<? extends X> callback) throws IOException {
    try (RandomAccessFile file = new RandomAccessFile(fastaFile, "r")) {
      try (FileChannel channel = file.getChannel()) {
        class Contig {

          private final int basesPerLine;
          private final ByteBuffer contig;
          private final String contigName;
          private final int length;

          Contig(ByteBuffer contig, int length, String contigName, int basesPerLine) {
            this.contig = contig;
            this.length = length;
            this.contigName = contigName;
            this.basesPerLine = basesPerLine;
          }

          private int accountForNewlines(int i) {
            return i + i / basesPerLine;
          }

          String get(int beginIndex, int endIndex) {
            StringBuilder builder = new StringBuilder();
            for (int end = Math.min(length, accountForNewlines(endIndex)),
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
            return length;
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
            new FastaFile() {

              @Override public int contigLength(String contig) {
                return Optional.ofNullable(chromosomes.get(contig))
                    .map(Contig::length)
                    .orElse(0);
              }

              @Override public Set<String> contigs() {
                return chromosomes.keySet();
              }

              @Override public String get(String contig, int beginIndex, int endIndex) {
                return Optional.ofNullable(chromosomes.get(contig))
                    .map(string -> string.get(beginIndex, endIndex))
                    .orElseThrow(() -> new IllegalArgumentException(
                        String.format("No such contig: \"%s\"", contig)));
              }
            });
      }
    }
  }
}

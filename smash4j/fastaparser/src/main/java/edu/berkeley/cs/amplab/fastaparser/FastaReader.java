package edu.berkeley.cs.amplab.fastaparser;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.Map;

public class FastaReader {

  public interface Callback<X> {

    public interface FastaFile {

      Map<String, Integer> info();

      Optional<String> get(String contig, int beginIndex, int endIndex);
    }

    X read(FastaFile fastaFile) throws Exception;
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
        return callback.read(
            new Callback.FastaFile() {

              class Contig {

                final CharSequence contig;
                final int basesPerLine;

                Contig(CharSequence contig, int basesPerLine) {
                  this.contig = contig;
                  this.basesPerLine = basesPerLine;
                }
              }

              private final Function<Contig, Integer> length =
                  new Function<Contig, Integer>() {
                    @Override public Integer apply(Contig contig) {
                      return contig.contig.length();
                    }
                  };

              private final Map<String, Contig> chromosomes = chromosomes();

              @Override
              public Map<String, Integer> info() {
                return Maps.transformValues(chromosomes, length);
              }

              @Override public Optional<String> get(
                  String contig, final int beginIndex, final int endIndex) {
                return Optional.fromNullable(chromosomes.get(contig))
                    .transform(
                        new Function<Contig, String>() {
                          @Override public String apply(Contig contig) {
                            CharSequence string = contig.contig.subSequence(
                                beginIndex + beginIndex / contig.basesPerLine,
                                endIndex + endIndex / contig.basesPerLine);
                            int length = string.length();
                            StringBuilder builder = new StringBuilder();
                            for (int i = 0; i < length; ++i) {
                              char c = string.charAt(i);
                              switch (c) {
                                case 'A': case 'C': case 'G': case 'T': case 'U': case 'R':
                                case 'Y': case 'K': case 'M': case 'S': case 'W': case 'B':
                                case 'D': case 'H': case 'V': case 'N': case 'X': case '-':
                                  builder.append(c);
                                  break;
                                case '\n':
                                  break;
                                default:
                                  throw new IllegalStateException(
                                      String.format("Illegal character: %c", c));
                              }
                            }
                            return builder.toString();
                          }
                        });
              }

              private Map<String, Contig> chromosomes() throws IOException {
                final ImmutableMap.Builder<String, Contig> chromosomes = ImmutableMap.builder();
                for (FastaIndex.Entry entry : index.entries()) {
                  long length = entry.length();
                  int bases = entry.bases();
                  chromosomes.put(
                      entry.name(),
                      new Contig(
                          StandardCharsets.UTF_8.decode(channel.map(
                              FileChannel.MapMode.READ_ONLY,
                              entry.offset(),
                              (length / bases) * entry.bytes() + length % bases)),
                          entry.bases()));
                }
                return chromosomes.build();
              }
            });
      }
    }
  }
}

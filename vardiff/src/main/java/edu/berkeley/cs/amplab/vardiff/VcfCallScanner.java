package edu.berkeley.cs.amplab.vardiff;

import com.google.common.base.Preconditions;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Iterables;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class VcfCallScanner implements CallScanner {

  private static class Indexer<X> implements Collector<X, Indexer<X>, Map<X, Integer>> {

    static <X> Indexer<X> create() {
      return new Indexer<>();
    }

    private final ImmutableMap.Builder<X, Integer> builder = ImmutableMap.builder();
    private int index = 0;

    private Indexer() {}

    @Override public BiConsumer<Indexer<X>, X> accumulator() {
      return (indexer, object) -> indexer.builder.put(object, index++);
    }

    @Override public Set<java.util.stream.Collector.Characteristics> characteristics() {
      return Collections.emptySet();
    }

    @Override public BinaryOperator<Indexer<X>> combiner() {
      return (lhs, rhs) -> {
        throw new UnsupportedOperationException();
      };
    }

    @Override public Function<Indexer<X>, Map<X, Integer>> finisher() {
      return indexer -> indexer.builder.build();
    }

    @Override public Supplier<Indexer<X>> supplier() {
      return Indexer::new;
    }
  }

  private static final Pattern
      HEADER_PATTERN = Pattern.compile(
          Stream.of("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
              .collect(Collectors.joining("\t"))),
      SAMPLE_PATTERN = Pattern.compile("\t([^\t]+)"),
      RECORD_PATTERN = Pattern.compile(StreamSupport
          .stream(
              Spliterators.spliteratorUnknownSize(
                  Iterables.cycle(Stream.of("([^\t]+)").collect(Collectors.toList())).iterator(),
                  Spliterator.IMMUTABLE),
              false)
          .limit(9)
          .collect(Collectors.joining("\t"))),
      ALT_PATTERN = Pattern.compile("([^,]+)(?:,|$)"),
      FORMAT_PATTERN = Pattern.compile("([^:]+)(?::|$)"),
      GENOTYPE_PATTERN = Pattern.compile("(\\p{Digit}+)([|/]|$)");

  public static VcfCallScanner create(File vcf) {
    return new VcfCallScanner(vcf, Optional.empty());
  }

  public static VcfCallScanner create(File vcf, String sampleId) {
    return new VcfCallScanner(vcf, Optional.of(sampleId));
  }

  private static Call call(String line, int sampleIndex) {
    Matcher matcher = RECORD_PATTERN.matcher(line);
    Preconditions.checkState(matcher.lookingAt());
    final String contig = matcher.group(1);
    final int position = Integer.parseInt(matcher.group(2));
    final String reference = matcher.group(4);
    List<String> alts = stream(ALT_PATTERN.matcher(matcher.group(5)))
        .map(result -> result.group(1))
        .collect(Collectors.toList());
    final List<String> alternates = Collections.singletonList(".").equals(alts)
        ? Collections.emptyList()
        : alts;
    Map<String, Integer> format = stream(FORMAT_PATTERN.matcher(matcher.group(9)))
        .map(result -> result.group(1))
        .collect(Indexer.create());
    List<String> call = stream(
            FORMAT_PATTERN.matcher(
                stream(matcher.usePattern(SAMPLE_PATTERN))
                    .map(result -> result.group(1))
                    .collect(Collectors.toList())
                    .get(sampleIndex)))
        .map(result -> result.group(1))
        .collect(Collectors.toList());
    boolean unphased = false, phased = false;
    final List<Integer> genotype = new ArrayList<>();
    for (
        Iterator<MatchResult>
            iterator = stream(GENOTYPE_PATTERN.matcher(call.get(format.get("GT")))).iterator();
        iterator.hasNext();) {
      MatchResult next = iterator.next();
      genotype.add(Integer.parseInt(next.group(1)));
      switch (next.group(2)) {
        case "|":
          phased = true;
          break;
        case "/":
          unphased = true;
          break;
        case "":
          break;
        default:
          throw new IllegalStateException();
      }
    }
    if (unphased && phased) {
      throw new IllegalStateException("Genotypes are either phased or unphased");
    }
    final Optional<Call.Phaseset> phaseset = phased
        ? Optional.of(
            Optional.ofNullable(format.get("PS"))
                .map(index -> Call.Phaseset.create(Integer.parseInt(call.get(index))))
                .orElse(Call.Phaseset.DEFAULT))
        : Optional.empty();
    return new Call() {

          @Override public List<String> alternates() {
            return alternates;
          }

          @Override public String contig() {
            return contig;
          }

          @Override public boolean equals(Object obj) {
            return HASH_CODE_AND_EQUALS.equals(this, obj);
          }

          @Override public List<Integer> genotype() {
            return genotype;
          }

          @Override public int hashCode() {
            return HASH_CODE_AND_EQUALS.hashCode(this);
          }

          @Override public Optional<Phaseset> phaseset() {
            return phaseset;
          }

          @Override public int position() {
            return position;
          }

          @Override public String reference() {
            return reference;
          }

          @Override public String toString() {
            return TO_STRING.apply(this);
          }
        };
  }

  private static Stream<MatchResult> stream(final Matcher matcher) {
    return StreamSupport.stream(
        Spliterators.spliteratorUnknownSize(
            new AbstractIterator<MatchResult>() {
              @Override protected MatchResult computeNext() {
                return matcher.find() ? matcher : endOfData();
              }
            },
            Spliterator.IMMUTABLE | Spliterator.NONNULL),
        false);
  }

  private final Optional<String> sampleId;
  private final File vcf;

  private VcfCallScanner(File vcf, Optional<String> sampleId) {
    this.vcf = vcf;
    this.sampleId = sampleId;
  }

  @Override
  public <X> X scan(Callback<? extends X> callback) throws IOException {
    try (BufferedReader in = new BufferedReader(new FileReader(vcf))) {
      List<String> header = new ArrayList<>();
      Stream<String> lines = in.lines().filter(line -> !(line.isEmpty() || line.startsWith("##")));
      Spliterator<String> spliterator = lines.spliterator();
      spliterator.tryAdvance(header::add);
      lines = StreamSupport.stream(spliterator, lines.isParallel());
      Matcher matcher = HEADER_PATTERN.matcher(Iterables.getOnlyElement(header));
      Preconditions.checkState(matcher.lookingAt());
      Map<String, Integer> index1 = stream(matcher.usePattern(SAMPLE_PATTERN))
          .map(result -> result.group(1))
          .collect(Indexer.create());
      int index = sampleId
          .map(id -> Optional.ofNullable(index1.get(id))
              .orElseThrow(() -> new IllegalStateException("Sample ID not present on header line")))
          .orElseGet(() -> {
            if (1 == index1.size()) {
              return 0;
            }
            throw new IllegalStateException("Sample ID required for multi-sample VCF file");
          });
      return callback.scan(lines.map(line -> call(line, index)));
    }
  }
}

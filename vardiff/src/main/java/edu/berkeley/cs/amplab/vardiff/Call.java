package edu.berkeley.cs.amplab.vardiff;

import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public interface Call extends Comparable<Call> {

  public static class Phaseset {

    public static final Phaseset DEFAULT = new Phaseset(Optional.empty());

    @SuppressWarnings("hiding")
    private static final HashCodeAndEquals<Phaseset> HASH_CODE_AND_EQUALS =
        HashCodeAndEquals.create(Phaseset.class, Phaseset::value);

    public static Phaseset create(int value) {
      return new Phaseset(Optional.of(value));
    }

    private final Optional<Integer> value;

    private Phaseset(Optional<Integer> value) {
      this.value = value;
    }

    @Override
    public boolean equals(Object obj) {
      return HASH_CODE_AND_EQUALS.equals(this, obj);
    }

    @Override
    public int hashCode() {
      return HASH_CODE_AND_EQUALS.hashCode(this);
    }

    @Override
    public String toString() {
      return value().map(Object::toString).orElse("DEFAULT");
    }

    public Optional<Integer> value() {
      return value;
    }
  }

  final Comparator<Call>
      COMPARATOR = Comparator.comparing(Call::contig).thenComparing(Call::position);

  final HashCodeAndEquals<Call> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
      Call.class,
      Call::alternates,
      Call::contig,
      Call::genotype,
      Call::phaseset,
      Call::position,
      Call::reference);

  final Function<Call, String> TO_STRING = call -> Stream
      .of(
          Stream.of("contig", call.contig()),
          Stream.of("position", call.position()),
          Stream.of("reference", call.reference()),
          Stream.of("alternates", call.alternates()),
          Stream.of("genotype", call.genotype()),
          Stream.of("phaseset", call.phaseset()))
      .map(stream -> stream.map(Object::toString).collect(Collectors.joining(" = ")))
      .collect(Collectors.joining(", "));

  List<String> alternates();

  @Override
  default int compareTo(Call rhs) {
    return COMPARATOR.compare(this, rhs);
  }

  String contig();

  default int end() {
    return position() + reference().length();
  }

  List<Integer> genotype();

  default boolean overlaps(Call rhs) {
    return Objects.equals(contig(), rhs.contig())
        && end() < rhs.position()
        && position() < rhs.end();
  }

  Optional<Phaseset> phaseset();

  int position();

  String reference();
}

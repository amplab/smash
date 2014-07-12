package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.ImmutableMap;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

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
      COMPARATOR = Comparator.comparing(Call::contig).thenComparing(Call::position),
      START_COMPARATOR = Comparator.comparing(Call::position),
      END_COMPARATOR = Comparator.comparing(call -> call.position() + call.reference().length());

  final HashCodeAndEquals<Call> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
      Call.class,
      Call::alternates,
      Call::contig,
      Call::genotype,
      Call::phaseset,
      Call::position,
      Call::reference);

  final Function<Call, String> TO_STRING = call -> ImmutableMap
      .builder()
      .put("contig", call.contig())
      .put("position", call.position())
      .put("reference", call.reference())
      .put("alternates", call.alternates())
      .put("genotype", call.genotype())
      .put("phaseset", call.phaseset())
      .build()
      .entrySet()
      .stream()
      .map(entry -> String.format("%s: %s", entry.getKey(), entry.getValue()))
      .collect(Collectors.joining(", "));

  List<String> alternates();

  @Override
  default int compareTo(Call rhs) {
    return COMPARATOR.compare(this, rhs);
  }

  String contig();

  List<Integer> genotype();

  Optional<Phaseset> phaseset();

  int position();

  String reference();
}

package edu.berkeley.cs.amplab.calldiff;

import com.google.api.services.genomics.model.Callset;
import com.google.common.collect.Iterables;

import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * An interface describing the information in a variant call used to determine if two sets of
 * calls are equivalent. There are two implementations: {@link VcfCallScanner} for using a sample
 * column from a VCF file as a set of calls, and {@link ApiCallScanner}, for using a {@link Callset}
 * as a set of calls.
 */
public interface Call {

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

  public enum Type {

    DELETION,
    HOM_REF,
    INSERTION,
    INVERSION,
    INDEL_OTHER,
    SNP;

    public static Type classify(Call call) {
      if (!call.genotype().stream().allMatch(Predicate.isEqual(0))) {
        List<String> alternates = call.alternates();
        if (1 == alternates.size()) {
          String
              reference = call.reference(),
              alternate = Iterables.getOnlyElement(alternates);
          int referenceSize = reference.length(),
              alternateSize = alternate.length();
          return 1 == referenceSize
              ? 1 == alternateSize
                  ? SNP
                  : INSERTION
              : 1 == alternateSize
                  ? DELETION
                  : Objects.equals(reference, new StringBuilder(alternate).reverse().toString())
                      ? INVERSION
                      : INDEL_OTHER;
        }
        return INDEL_OTHER;
      }
      return HOM_REF;
    }
  }

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

  String contig();

  default int end() {
    return position() + reference().length();
  }

  List<Integer> genotype();

  default boolean overlaps(Call rhs) {
    return Objects.equals(contig(), rhs.contig())
        && position() < rhs.end()
        && rhs.position() < end();
  }

  Optional<Phaseset> phaseset();

  int position();

  String reference();

  default Type type() {
    return Type.classify(this);
  }
}

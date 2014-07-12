package edu.berkeley.cs.amplab.vardiff;

import static org.junit.Assert.assertTrue;

import com.google.common.base.Throwables;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Lists;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.Random;
import java.util.function.BiPredicate;
import java.util.function.Supplier;

public class TestCall implements Call {

  public static void checkNonOverlapping(List<Call> calls) {
    checkNonOverlapping(calls, calls, (x, y) -> x == y);
  }

  public static void checkNonOverlapping(List<Call> lhs, List<Call> rhs) {
    checkNonOverlapping(lhs, rhs, (x, y) -> false);
  }

  private static void checkNonOverlapping(List<Call> lhs, List<Call> rhs,
      BiPredicate<Call, Call> predicate) {
    for (Call call1 : lhs) {
      for (Call call2 : rhs) {
        assertTrue(predicate.test(call1, call2) || !overlaps(call1, call2));
      }
    }
  }

  public static TestCall create(
      String contig,
      int position,
      String reference,
      List<String> alternates,
      List<Integer> genotype) {
    return create(
        contig,
        position,
        reference,
        alternates,
        genotype,
        Optional.empty());
  }

  public static TestCall create(
      String contig,
      int position,
      String reference,
      List<String> alternates,
      List<Integer> genotype,
      Optional<Phaseset> phaseset) {
    return new TestCall(
        contig,
        position,
        reference,
        alternates,
        genotype,
        phaseset);
  }

  public static TestCall create(
      String contig,
      int position,
      String reference,
      List<String> alternates,
      List<Integer> genotype,
      Phaseset phaseset) {
    return create(
        contig,
        position,
        reference,
        alternates,
        genotype,
        Optional.of(phaseset));
  }

  public static boolean overlaps(Call call1, Call call2) {
    int start1 = call1.position(),
        start2 = call2.position();
    return Objects.equals(call1.contig(), call2.contig())
        && start1 < start2 + call2.reference().length()
        && start2 < start1 + call1.reference().length();
  }

  public static List<Call> randomCalls(
      Random random,
      String contig,
      int contigLength,
      int maxCallLength,
      int numberOfCalls) {
    return randomCalls(
        random,
        contig,
        contigLength,
        maxCallLength,
        numberOfCalls,
        () -> Optional.empty());
  }

  public static List<Call> randomCalls(
      Random random,
      String contig,
      int contigLength,
      int maxCallLength,
      int numberOfCalls,
      Supplier<Optional<Phaseset>> phaseset) {
    try {
      ImmutableSortedSet.Builder<Call> calls = ImmutableSortedSet.orderedBy(Call.START_COMPARATOR);
      for (int i = 0; i < numberOfCalls; ++i) {
        int length = 1 == maxCallLength ? 1 : random.nextInt(maxCallLength - 1) + 1,
            max = contigLength - length,
            start = 0 == max ? 1 : random.nextInt(max);
        calls.add(create(
            contig,
            start,
            TestReference.reader().read(
                reference -> reference.get(contig, start - 1, start + length - 1)),
            Collections.emptyList(),
            Arrays.asList(0, 0),
            phaseset.get()));
      }
      return Lists.newArrayList(calls.build());
    } catch (IOException e) {
      throw Throwables.propagate(e);
    }
  }

  private final List<String> alternates;
  private final String contig;
  private final List<Integer> genotype;
  private final Optional<Phaseset> phaseset;
  private final int position;
  private final String reference;

  private TestCall(
      String contig,
      int position,
      String reference,
      List<String> alternates,
      List<Integer> genotype,
      Optional<Phaseset> phaseset) {
    this.contig = contig;
    this.position = position;
    this.reference = reference;
    this.alternates = alternates;
    this.genotype = genotype;
    this.phaseset = phaseset;
  }

  @Override
  public List<String> alternates() {
    return alternates;
  }

  @Override
  public String contig() {
    return contig;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public List<Integer> genotype() {
    return genotype;
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  @Override
  public Optional<Phaseset> phaseset() {
    return phaseset;
  }

  @Override
  public int position() {
    return position;
  }

  @Override
  public String reference() {
    return reference;
  }

  @Override
  public String toString() {
    return TO_STRING.apply(this);
  }
}

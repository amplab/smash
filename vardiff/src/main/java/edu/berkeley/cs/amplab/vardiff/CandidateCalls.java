package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.Sets;

import com.beust.jcommander.internal.Lists;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

public class CandidateCalls {

  private static final HashCodeAndEquals<CandidateCalls>
      HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
          CandidateCalls.class,
          CandidateCalls::lhs,
          CandidateCalls::rhs);

  public static List<CandidateCalls> createCandidates(List<Call> lhsCalls, List<Call> rhsCalls) {
    List<List<Call>>
        lhsSubsets = allNonOverlappingSubsets(lhsCalls),
        rhsSubsets = allNonOverlappingSubsets(rhsCalls);
    List<CandidateCalls> candidates = new ArrayList<>();
    for (List<Call> lhs : lhsSubsets) {
      for (List<Call> rhs : rhsSubsets) {
        candidates.add(new CandidateCalls(lhs, rhs));
      }
    }
    Collections.sort(
        candidates,
        (lhs, rhs) -> Integer.compare(
            rhs.lhs().size() + rhs.rhs().size(),
            lhs.lhs().size() + lhs.rhs().size()));
    return candidates;
  }

  static List<List<Call>> allNonOverlappingSubsets(List<Call> calls) {
    Set<Call> sortedByEnd = Sets.newTreeSet(Call.END_COMPARATOR);
    Stream.Builder<Set<Call>> streamBuilder = Stream.builder();
    for (sortedByEnd.addAll(calls); !sortedByEnd.isEmpty();) {
      Set<Call>
          selected = Sets.newTreeSet(Call.START_COMPARATOR),
          copy = Sets.newTreeSet(Call.END_COMPARATOR);
      for (copy.addAll(sortedByEnd); !copy.isEmpty();) {
        Iterator<Call> iterator = copy.iterator();
        Call call = iterator.next();
        int start = call.position(),
            end = end(call);
        iterator.remove();
        for (selected.add(call); iterator.hasNext();) {
          Call next = iterator.next();
          if (start < end(next) && next.position() < end) {
            iterator.remove();
          }
        }
      }
      streamBuilder.add(selected);
      sortedByEnd.removeAll(selected);
    }
    return Stream
        .concat(
            streamBuilder.build()
                .flatMap(set -> Sets.powerSet(set).stream())
                .filter(set -> !set.isEmpty())
                .map(Lists::newArrayList),
            Stream.of(Collections.<Call>emptyList()))
        .collect(
            ArrayList::new,
            (list, object) -> list.add(object),
            (lhs, rhs) -> { throw new UnsupportedOperationException(); });
  }

  private static int end(Call call) {
    int start = call.position();
    return start + call.reference().length();
  }

  private final List<Call> lhs;
  private final List<Call> rhs;

  private CandidateCalls(List<Call> lhs, List<Call> rhs) {
    this.lhs = lhs;
    this.rhs = rhs;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public List<Call> lhs() {
    return lhs;
  }

  public List<Call> rhs() {
    return rhs;
  }

  public int size() {
    return lhs().size() + rhs.size();
  }
}

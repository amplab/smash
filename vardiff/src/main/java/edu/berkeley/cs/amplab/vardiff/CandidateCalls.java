package edu.berkeley.cs.amplab.vardiff;

import com.google.common.collect.Sets;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CandidateCalls {

  private static final HashCodeAndEquals<CandidateCalls> HASH_CODE_AND_EQUALS =
      HashCodeAndEquals.create(CandidateCalls.class,
          CandidateCalls::contig,
          CandidateCalls::start,
          CandidateCalls::end,
          CandidateCalls::lhs,
          CandidateCalls::rhs);

  public static CandidateCalls create(
      String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
    return new CandidateCalls(contig, start, end, lhs, rhs);
  }

  public static Stream<CandidateCalls> createCandidates(Window window) {
    List<Call>
        lhsWindow = window.lhs(),
        rhsWindow = window.rhs();
    Comparator<Call>
        lhsComparator = Comparator.comparing(lhsWindow.stream().collect(Indexer.create())::get),
        rhsComparator = Comparator.comparing(rhsWindow.stream().collect(Indexer.create())::get);
    List<Set<Call>>
        lhsPowerSet = nonOverlappingSubsets(lhsWindow),
        rhsPowerSet = nonOverlappingSubsets(rhsWindow);
    List<CandidateCalls>
        candidates = new ArrayList<>(lhsPowerSet.size() * rhsPowerSet.size());
    for (Set<Call> lhsSet : lhsPowerSet) {
      List<Call> lhs = new ArrayList<>(lhsSet.size());
      lhs.addAll(lhsSet);
      Collections.sort(lhs, lhsComparator);
      for (Set<Call> rhsSet : rhsPowerSet) {
        List<Call> rhs = new ArrayList<>(rhsSet.size());
        rhs.addAll(rhsSet);
        Collections.sort(rhs, rhsComparator);
        candidates.add(create(window.contig(), window.start(), window.end(), lhs, rhs));
      }
    }
    Collections.sort(candidates, Comparator.comparing(CandidateCalls::size).reversed());
    return candidates.stream();
  }

  private static List<Set<Call>> nonOverlappingSubsets(List<Call> list) {
    return Sets.powerSet(Sets.newHashSet(list))
        .stream()
        .filter(
            calls -> {
              for (Call lhs : calls) {
                for (Call rhs : calls) {
                  if (lhs != rhs && lhs.overlaps(rhs)) {
                    return false;
                  }
                }
              }
              return true;
            })
        .collect(Collectors.toList());
  }

  private final String contig;
  private final List<Call> lhs, rhs;
  private final int start, end;

  private CandidateCalls(String contig, int start, int end, List<Call> lhs, List<Call> rhs) {
    this.contig = contig;
    this.start = start;
    this.end = end;
    this.lhs = lhs;
    this.rhs = rhs;
  }

  public String contig() {
    return contig;
  }

  public int end() {
    return end;
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

  int size() {
    return lhs().size() + rhs().size();
  }

  public int start() {
    return start;
  }

  public boolean generatesSameSetOfHaplotypes(FastaReader.FastaFile reference) {
    String contig = contig();
    int start = start(), end = end();
    return Objects.equals(
        HaplotypeGenerator.generateHaplotypes(reference, contig, lhs(), start, end),
        HaplotypeGenerator.generateHaplotypes(reference, contig, rhs(), start, end));
  }
}

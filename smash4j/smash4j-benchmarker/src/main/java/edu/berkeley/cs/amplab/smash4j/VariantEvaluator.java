package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Objects;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Iterables;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;

import com.beust.jcommander.internal.Lists;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VariantEvaluator {

  public static class Builder {

    private int maxIndelLength = 50;
    private int maxSvBreakpointDistance = 100;
    private int maxVariantLengthDifference = 100;
    private Optional<FastaReader.Callback.FastaFile> reference = Optional.absent();
    private int rescueWindowSize = 50;

    private Builder() {}

    public VariantEvaluator build() {
      return new VariantEvaluator(
          reference,
          maxIndelLength,
          maxSvBreakpointDistance,
          maxVariantLengthDifference,
          rescueWindowSize);
    }

    public Builder setMaxIndelLength(int maxIndelLength) {
      this.maxIndelLength = maxIndelLength;
      return this;
    }

    public Builder setMaxSvBreakpointDistance(int maxSvBreakpointDistance) {
      this.maxSvBreakpointDistance = maxSvBreakpointDistance;
      return this;
    }

    public Builder setMaxVariantLengthDifference(int maxVariantLengthDifference) {
      this.maxVariantLengthDifference = maxVariantLengthDifference;
      return this;
    }

    public Builder setReference(FastaReader.Callback.FastaFile reference) {
      this.reference = Optional.of(reference);
      return this;
    }

    public Builder setRescueWindowSize(int rescueWindowSize) {
      this.rescueWindowSize = rescueWindowSize;
      return this;
    }
  }

  public static class ContigStats {

    static class Builder {

      private final NavigableMap<Long, VariantProto> true_var;
      private final NavigableMap<Long, VariantProto> pred_var;
      private final ImmutableList.Builder<Long> intersect_good;
      private final NavigableSet<Long> false_positives;
      private final NavigableSet<Long> false_negatives;
      private final GenotypeConcordance.Builder genotype_concordance;
      private Multimap<VariantType, Long> intersect_bad;
      private Optional<Multiset<VariantType>> known_fp = Optional.absent();
      private Optional<Multiset<VariantType>> calls_at_known_fp = Optional.absent();
      private Optional<List<Long>> known_fp_calls_positions = Optional.absent();

      private Builder(
          NavigableMap<Long, VariantProto> true_var,
          NavigableMap<Long, VariantProto> pred_var,
          ImmutableList.Builder<Long> intersect_good,
          NavigableSet<Long> false_positives,
          NavigableSet<Long> false_negatives,
          GenotypeConcordance.Builder genotype_concordance) {
        this.true_var = true_var;
        this.pred_var = pred_var;
        this.intersect_good = intersect_good;
        this.false_positives = false_positives;
        this.false_negatives = false_negatives;
        this.genotype_concordance = genotype_concordance;
      }

      Builder rectify(FastaReader.Callback.FastaFile reference, int rescueWindowSize) {
        Collection<Long> locs_to_rescue = Lists.newArrayList(false_negatives);
        for (Long loc : locs_to_rescue) {
          if (false_negatives.contains(loc)) {
//          new_tp,rm_fp,rescued_vars = rescue_mission(self.false_negatives,self.false_positives,self.true_positives,loc,ref,window)
//          for t in VARIANT_TYPE:
//              # seemingly odd accounting. The number of predicted variants *changes* as a result of rescuing.
//              # e.g. 2 predicted FPs are in fact 1 FN. So
//              #  -- remove 2 predicted variants
//              #  -- remove 2 false positives
//              #  -- remove 1 false negative
//              #  -- add 1 true positive
//              self.num_pred[t] -= rm_fp[t]
//              self.num_fp[t] -= rm_fp[t]
//              self.num_pred[t] += new_tp[t]
//              self.num_fn[t] -= new_tp[t]
//              self.num_tp[t] += new_tp[t]
//          for v in rescued_vars:
//              self.rescued_vars._add_variant(v)
          }
        }
        return this;
      }

      Builder setIntersectBad(Multimap<VariantType, Long> intersect_bad) {
        this.intersect_bad = intersect_bad;
        return this;
      }

      Builder setKnownFp(Multiset<VariantType> known_fp) {
        this.known_fp = Optional.of(known_fp);
        return this;
      }

      Builder setCallsAtKnownFp(Multiset<VariantType> calls_at_known_fp) {
        this.calls_at_known_fp = Optional.of(calls_at_known_fp);
        return this;
      }

      Builder setKnownFpCallsPositions(List<Long> known_fp_calls_positions) {
        this.known_fp_calls_positions = Optional.of(known_fp_calls_positions);
        return this;
      }

      ContigStats build() {
        return new ContigStats();
      }
    }

    static Builder builder(
        NavigableMap<Long, VariantProto> true_var,
        NavigableMap<Long, VariantProto> pred_var,
        ImmutableList.Builder<Long> intersect_good,
        NavigableSet<Long> false_positives,
        NavigableSet<Long> false_negatives,
        GenotypeConcordance.Builder genotype_concordance) {
      return new Builder(
          true_var,
          pred_var,
          intersect_good,
          false_positives,
          false_negatives,
          genotype_concordance);
    }

    private ContigStats() {}
  }

  public enum Genotype {

    HET,
    HOM_REF,
    HOM_VAR,
    NO_CALL;

    private static final Function<String, Genotype> PARSE_GENOTYPE =
        new Function<String, Genotype>() {

          private final Pattern splitPattern = Pattern.compile("(0|(?:[1-9][0-9]*?))(?:[/|]|$)");
          private final Integer zero = 0;

          @Override public Genotype apply(String genotype) {
            Integer i = null;
            for (Matcher matcher = splitPattern.matcher(genotype); matcher.matches();) {
              Integer j = Integer.parseInt(matcher.group(1));
              if (null == i) {
                i = j;
              } else if (!i.equals(j)) {
                return HET;
              }
            }
            return null == i ? NO_CALL : zero.equals(i) ? HOM_REF : HOM_VAR;
          }
        };

    public static Genotype getGenotype(VariantProto variant) {
      return GenotypeExtractor.INSTANCE.getGenotype(variant).transform(PARSE_GENOTYPE).or(NO_CALL);
    }
  }

  public static class GenotypeConcordance {

    static class Builder {

      private static final
          Function<
              Map<Genotype, ImmutableMultiset.Builder<Genotype>>,
              Map<Genotype, Multiset<Genotype>>>
          BUILD =
          new Function<
              Map<Genotype, ImmutableMultiset.Builder<Genotype>>,
              Map<Genotype, Multiset<Genotype>>>() {

            private final Function<ImmutableMultiset.Builder<Genotype>, Multiset<Genotype>> build =
                new Function<ImmutableMultiset.Builder<Genotype>, Multiset<Genotype>>() {
                  @Override public Multiset<Genotype> apply(
                      ImmutableMultiset.Builder<Genotype> builder) {
                    return builder.build();
                  }
                };

            @Override public Map<Genotype, Multiset<Genotype>> apply(
                Map<Genotype, ImmutableMultiset.Builder<Genotype>> map) {
              return Maps.transformValues(map, build);
            }
          };

      private final
          Map<VariantType, Map<Genotype, ImmutableMultiset.Builder<Genotype>>> concordance;

      private Builder() {
        ImmutableMap.Builder<VariantType, Map<Genotype, ImmutableMultiset.Builder<Genotype>>>
            concordance = ImmutableMap.builder();
        for (VariantType variantType : VariantType.values()) {
          ImmutableMap.Builder<Genotype, ImmutableMultiset.Builder<Genotype>>
              builder = ImmutableMap.builder();
          for (Genotype genotype : Genotype.values()) {
            builder.put(genotype, ImmutableMultiset.<Genotype>builder());
          }
          concordance.put(variantType, builder.build());
        }
        this.concordance = concordance.build();
      }

      GenotypeConcordance build() {
        return new GenotypeConcordance(Maps.transformValues(concordance, BUILD));
      }

      Builder increment(
          VariantType variantType,
          Genotype trueGenotype,
          Genotype predictedGenotype) {
        concordance.get(variantType).get(trueGenotype).add(predictedGenotype);
        return this;
      }
    }

    static Builder builder() {
      return new Builder();
    }

    private final Map<VariantType, Map<Genotype, Multiset<Genotype>>> concordance;

    private GenotypeConcordance(Map<VariantType, Map<Genotype, Multiset<Genotype>>> concordance) {
      this.concordance = concordance;
    }

    public int get(VariantType variantType, Genotype trueGenotype, Genotype predictedGenotype) {
      return concordance.get(variantType).get(trueGenotype).count(predictedGenotype);
    }
  }

  public enum VariantType {

    INDEL_DELETION(false, true, false, false, true, false, false),
    INDEL_INSERTION(false, true, false, true, false, false, false),
    INDEL_INVERSION(false, true, false, false, false, true, false),
    INDEL_OTHER(false, true, false, false, false, false, true),
    SNP(true, false, false, false, false, false, false),
    SV_DELETION(false, false, true, false, true, false, false),
    SV_INSERTION(false, false, true, true, false, false, false),
    SV_OTHER(false, false, true, false, false, false, true);

    private static String getFirstAlt(VariantProto variant) {
      return Iterables.getOnlyElement(variant.getAlternateBasesList());
    }

    public static VariantType getType(VariantProto variant) {
      return isSnp(variant)
          ? SNP
          : isSv(variant)
              ? hasSingleAlt(variant)
                  ? isInsertion(variant)
                      ? SV_INSERTION
                      : isDeletion(variant) ? SV_DELETION : SV_OTHER
                  : SV_OTHER
              : hasSingleAlt(variant)
                  ? isInsertion(variant)
                      ? INDEL_INSERTION
                      : isDeletion(variant)
                          ? INDEL_DELETION
                          : isInversion(variant) ? INDEL_INVERSION : INDEL_OTHER
                  : INDEL_OTHER;
    }

    private static boolean hasSingleAlt(VariantProto variant) {
      return 1 == variant.getAlternateBasesCount();
    }

    private static boolean isDeletion(VariantProto variant) {
      return getFirstAlt(variant).length() < variant.getReferenceBases().length();
    }

    private static boolean isInsertion(VariantProto variant) {
      return variant.getReferenceBases().length() < getFirstAlt(variant).length();
    }

    private static boolean isInversion(VariantProto variant) {
      String ref = variant.getReferenceBases(), alt = getFirstAlt(variant);
      int length = ref.length();
      if (length == alt.length()) {
        for (int i = 0; i < length; ++i) {
          if (ref.charAt(i) != alt.charAt(length - i - 1)) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    private static boolean isSnp(VariantProto variant) {
      if (1 == variant.getReferenceBases().length()) {
        for (String alternateBases : variant.getAlternateBasesList()) {
          if (1 != alternateBases.length()) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    private static boolean isSv(VariantProto variant) {
      for (VariantProto.Multimap.Entry entry : variant.getInfo().getEntryList()) {
        if ("SVTYPE".equals(entry.getKey())) {
          return true;
        }
      }
      return false;
    }

    private final boolean isDeletion;
    private final boolean isIndel;
    private final boolean isInsertion;
    private final boolean isInversion;
    private final boolean isOther;
    private final boolean isSnp;
    private final boolean isSv;

    private VariantType(
        boolean isSnp,
        boolean isIndel,
        boolean isSv,
        boolean isInsertion,
        boolean isDeletion,
        boolean isInversion,
        boolean isOther) {
      this.isSnp = isSnp;
      this.isIndel = isIndel;
      this.isSv = isSv;
      this.isInsertion = isInsertion;
      this.isDeletion = isDeletion;
      this.isInversion = isInversion;
      this.isOther = isOther;
    }

    public boolean isDeletion() {
      return isDeletion;
    }

    public boolean isIndel() {
      return isIndel;
    }

    public boolean isInsertion() {
      return isInsertion;
    }

    public boolean isInversion() {
      return isInversion;
    }

    public boolean isOther() {
      return isOther;
    }

    public boolean isSnp() {
      return isSnp;
    }

    public boolean isSv() {
      return isSv;
    }
  }

  private static final Function<VariantProto, String> GET_CONTIG =
      new Function<VariantProto, String>() {
        @Override public String apply(VariantProto variant) {
          return variant.getContig();
        }
      };

  private static final Function<VariantProto, Long> GET_POSITION =
      new Function<VariantProto, Long>() {
        @Override public Long apply(VariantProto variant) {
          return variant.getPosition();
        }
      };

  public static Builder builder() {
    return new Builder();
  }

  private static <X extends Comparable<? super X>, Y> NavigableMap<X, Y> uniqueIndex(
      Iterable<Y> iterable, Function<? super Y, X> function) {
    NavigableMap<X, Y> result = new TreeMap<>();
    for (Y object : iterable) {
      result.put(function.apply(object), object);
    }
    return result;
  }

  private final int maxIndelLength;
  private final int maxSvBreakpointDistance;
  private final int maxVariantLengthDifference;
  private final Optional<FastaReader.Callback.FastaFile> reference;
  private final int rescueWindowSize;

  private VariantEvaluator(
      Optional<FastaReader.Callback.FastaFile> reference,
      int maxIndelLength,
      int maxSvBreakpointDistance,
      int maxVariantLengthDifference,
      int rescueWindowSize) {
    this.reference = reference;
    this.maxIndelLength = maxIndelLength;
    this.maxSvBreakpointDistance = maxSvBreakpointDistance;
    this.maxVariantLengthDifference = maxVariantLengthDifference;
    this.rescueWindowSize = rescueWindowSize;
  }

  public Map<String, ContigStats> evaluate(
      Iterable<VariantProto> trueVariants,
      Iterable<VariantProto> predictedVariants,
      Iterable<VariantProto> knownFalsePositives) {
    ListMultimap<String, VariantProto>
        trueVars = Multimaps.index(trueVariants, GET_CONTIG),
        predVars = Multimaps.index(predictedVariants, GET_CONTIG),
        knownFps = Multimaps.index(knownFalsePositives, GET_CONTIG);
    ImmutableMap.Builder<String, ContigStats> result = ImmutableMap.builder();
    for (String contig :
        Sets.union(Sets.union(trueVars.keySet(), predVars.keySet()), knownFps.keySet())) {
      result.put(
          contig,
          evaluate(
              contig,
              uniqueIndex(trueVars.get(contig), GET_POSITION),
              uniqueIndex(predVars.get(contig), GET_POSITION),
              uniqueIndex(knownFps.get(contig), GET_POSITION)));
    }
    return result.build();
  }

  @SafeVarargs
  private static <X extends Comparable<? super X>> NavigableSet<X> intersection(
      Set<? extends X>... sets) {
    NavigableSet<X> intersection = new TreeSet<>();
    Collection<Set<? extends X>> collection = Arrays.asList(sets);
    for (Set<? extends X> set : Iterables.limit(collection, 1)) {
      intersection.addAll(set);
    }
    for (Set<? extends X> set : Iterables.skip(collection, 1)) {
      intersection.retainAll(set);
    }
    return intersection;
  }

  private static <X extends Comparable<? super X>> NavigableSet<X> difference(
      Set<? extends X> lhs, Set<? extends X> rhs) {
    NavigableSet<X> difference = new TreeSet<>();
    difference.addAll(lhs);
    difference.removeAll(rhs);
    return difference;
  }

  private ContigStats evaluate(
      String contig,
      NavigableMap<Long, VariantProto> true_var,
      NavigableMap<Long, VariantProto> pred_var,
      NavigableMap<Long, VariantProto> known_fp) {
    NavigableSet<Long>
        true_loc = true_var.navigableKeySet(),
        pred_loc = pred_var.navigableKeySet(),
        false_positives = difference(pred_loc, true_loc),
        false_negatives = difference(true_loc, pred_loc);
    GenotypeConcordance.Builder
        genotype_concordance = GenotypeConcordance.builder();
    ImmutableList.Builder<Long>
        intersect_good = ImmutableList.builder(),
        intersect_bad = ImmutableList.builder(),
        known_fp_calls_positions = ImmutableList.builder();
    ImmutableMultimap.Builder<VariantType, Long>
        intersect_bad_dict = ImmutableMultimap.builder();
    ImmutableMultiset.Builder<VariantType>
        calls_at_known_fp = ImmutableMultiset.builder(),
        all_known_fp = ImmutableMultiset.builder();
    boolean
        hasKnownFp = !known_fp.isEmpty();
    for (Long loc : intersection(pred_loc, true_loc)) {
      VariantProto
          t = true_var.get(loc),
          p = pred_var.get(loc);
      VariantType
          vartype = VariantType.getType(t);
      if (vartype == VariantType.getType(p)
          && (vartype.isSv() || t.getAlternateBasesList().equals(p.getAlternateBasesList()))) {
        intersect_good.add(loc);
        genotype_concordance.increment(
            vartype,
            Genotype.getGenotype(t),
            Genotype.getGenotype(p));
      } else {
        intersect_bad.add(loc);
        intersect_bad_dict.put(vartype, loc);
      }
    }
    if (hasKnownFp) {
      for (Long loc : intersection(pred_loc, known_fp.navigableKeySet())) {
        VariantProto
            k = known_fp.get(loc),
            p = pred_var.get(loc);
        VariantType
            vartype = VariantType.getType(k);
        if (Objects.equal(k.getReferenceBases(), p.getReferenceBases())) {
          calls_at_known_fp.add(vartype);
          known_fp_calls_positions.add(loc);
        }
        all_known_fp.add(vartype);
      }
    }
    for (Long loc : false_negatives) {
      VariantProto t = true_var.get(loc);
      VariantType vartype = VariantType.getType(t);
      if (vartype.isSv()) {
        Optional<Long> optional = structuralMatch(t, pred_var);
        if (optional.isPresent()) {
          Long match = optional.get();
          if (false_positives.contains(match)) {
            intersect_good.add(loc);
            false_positives.remove(match);
            false_negatives.remove(loc);
            genotype_concordance.increment(
                vartype,
                Genotype.getGenotype(t),
                Genotype.getGenotype(pred_var.get(match)));
          }
        }
      }
    }
    Collection<Long> ib = intersect_bad.build();
    false_positives.addAll(ib);
    false_negatives.addAll(ib);
    ContigStats.Builder variant_stats = ContigStats.builder(
        true_var,
        pred_var,
        intersect_good,
        false_positives,
        false_negatives,
        genotype_concordance);
    if (reference.isPresent()) {
      variant_stats.rectify(reference.get(), rescueWindowSize);
    }
    if (hasKnownFp) {
      variant_stats
          .setKnownFp(all_known_fp.build())
          .setCallsAtKnownFp(calls_at_known_fp.build())
          .setKnownFpCallsPositions(known_fp_calls_positions.build());
    }
    return variant_stats
        .setIntersectBad(intersect_bad_dict.build())
        .build();
  }

  private Optional<Long> structuralMatch(
      final VariantProto trueVariant,
      NavigableMap<Long, VariantProto> predictedVariants) {
    final long position = trueVariant.getPosition();
    Collection<VariantProto> candidates = FluentIterable
        .from(predictedVariants
            .subMap(
                position - maxSvBreakpointDistance, true,
                position + maxSvBreakpointDistance, true)
            .values())
        .filter(
            new Predicate<VariantProto>() {

              private final VariantType type = VariantType.getType(trueVariant);

              @Override public boolean apply(VariantProto predictedVariant) {
                if (type == VariantType.getType(predictedVariant)) {
                  int trueReferenceLength = trueVariant.getReferenceBases().length(),
                      predictedReferenceLength = predictedVariant.getReferenceBases().length();
                  for (String trueAlt : trueVariant.getAlternateBasesList()) {
                    int trueGain = trueReferenceLength - trueAlt.length();
                    for (String predictedAlt : predictedVariant.getAlternateBasesList()) {
                      int predGain = predictedReferenceLength - predictedAlt.length();
                      if (trueGain - maxVariantLengthDifference < predGain
                          && predGain < trueGain + maxVariantLengthDifference) {
                        return true;
                      }
                    }
                  }
                }
                return false;
              }
            })
        .toList();
    Queue<VariantProto> priorityQueue = new PriorityQueue<>(Math.max(candidates.size(), 1),
        new Comparator<VariantProto>() {
          @Override public int compare(VariantProto lhs, VariantProto rhs) {
            return Long.compare(
                Math.abs(position - lhs.getPosition()),
                Math.abs(position - rhs.getPosition()));
          }
        });
    priorityQueue.addAll(candidates);
    return Optional.fromNullable(priorityQueue.peek()).transform(GET_POSITION);
  }
}

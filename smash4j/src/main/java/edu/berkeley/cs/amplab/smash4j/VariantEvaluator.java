package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.repackaged.com.google.common.base.Joiner;
import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Ordering;

import edu.berkeley.cs.amplab.smash4j.fasta.FastaReader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class VariantEvaluator {

  public static class Builder {

    private int maxIndelSize = 50;
    private int maxSvBreakpointDistance = 100;
    private int maxVariantLengthDifference = 100;
    private Optional<FastaReader.FastaFile> reference = Optional.absent();
    private int rescueWindowSize = 50;

    private Builder() {}

    public VariantEvaluator build() {
      return new VariantEvaluator(
          reference,
          maxSvBreakpointDistance,
          maxVariantLengthDifference,
          rescueWindowSize,
          VariantType.getType(maxIndelSize));
    }

    public Builder setMaxIndexSize(int maxIndelSize) {
      this.maxIndelSize = maxIndelSize;
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

    public Builder setReference(FastaReader.FastaFile reference) {
      this.reference = Optional.of(reference);
      return this;
    }

    public Builder setRescueWindowSize(int rescueWindowSize) {
      this.rescueWindowSize = rescueWindowSize;
      return this;
    }
  }

  public static class ContigStats {

    static ContigStats create(
        String contig,
        NavigableMap<Integer, Variant> trueVariants,
        NavigableMap<Integer, Variant> predictedVariants,
        Iterable<Integer> truePositiveLocations,
        Iterable<Integer> falsePositiveLocations,
        Iterable<Integer> falseNegativeLocations,
        Iterable<Integer> incorrectPredictions,
        GenotypeConcordance concordance,
        Function<Variant, VariantType> getType) {
      NavigableMap<Integer, Variant>
          truePositives = filter(predictedVariants, truePositiveLocations),
          falsePositives = filter(predictedVariants, falsePositiveLocations),
          falseNegatives = filter(trueVariants, falseNegativeLocations);
      Multimap<VariantType, Integer> incorrectPredictionsMultimap = ArrayListMultimap.create();
      for (Integer location : incorrectPredictions) {
        incorrectPredictionsMultimap.put(getType.apply(predictedVariants.get(location)), location);
      }
      return new ContigStats(
          contig,
          trueVariants,
          predictedVariants,
          truePositives,
          falsePositives,
          falseNegatives,
          incorrectPredictionsMultimap,
          concordance,
          getType);
    }

    private static <X extends Comparable<? super X>, Y> NavigableMap<X, Y>
        filter(Map<? super X, ? extends Y> map, Iterable<? extends X> keys) {
      NavigableMap<X, Y> filtered = new TreeMap<>();
      for (X key : keys) {
        filtered.put(key, map.get(key));
      }
      return filtered;
    }

    private static String toString(Map<Integer, ?> map) {
      return Joiner.on(", ").join(map.keySet());
    }

    private static String toString(Optional<? extends Map<Integer, ?>> optional) {
      return optional.isPresent() ? toString(optional.get()) : "absent";
    }

    private Optional<NavigableMap<Integer, Variant>>
        allKnownFalsePositives = Optional.absent(),
        correctKnownFalsePositives = Optional.absent(),
        knownFalsePositives = Optional.absent(),
        rescuedVariants = Optional.absent();
    private final GenotypeConcordance concordance;
    private final String contig;
    private final NavigableMap<Integer, Variant> falseNegatives;
    private final NavigableMap<Integer, Variant> falsePositives;
    private final Multimap<VariantType, Integer> incorrectPredictions;
    private final Function<Variant, VariantType> getType;
    private final NavigableMap<Integer, Variant> predictedVariants;
    private final NavigableMap<Integer, Variant> truePositives;
    private final NavigableMap<Integer, Variant> trueVariants;

    private ContigStats(
        String contig,
        NavigableMap<Integer, Variant> trueVariants,
        NavigableMap<Integer, Variant> predictedVariants,
        NavigableMap<Integer, Variant> truePositives,
        NavigableMap<Integer, Variant> falsePositives,
        NavigableMap<Integer, Variant> falseNegatives,
        Multimap<VariantType, Integer> incorrectPredictions,
        GenotypeConcordance concordance,
        Function<Variant, VariantType> getType) {
      this.contig = contig;
      this.trueVariants = trueVariants;
      this.predictedVariants = predictedVariants;
      this.truePositives = truePositives;
      this.falsePositives = falsePositives;
      this.falseNegatives = falseNegatives;
      this.incorrectPredictions = incorrectPredictions;
      this.concordance = concordance;
      this.getType = getType;
    }

    public Optional<NavigableMap<Integer, Variant>> allKnownFalsePositives() {
      return allKnownFalsePositives;
    }

    public GenotypeConcordance concordance() {
      return concordance;
    }

    public String contig() {
      return contig;
    }

    public Optional<NavigableMap<Integer, Variant>> correctKnownFalsePositives() {
      return correctKnownFalsePositives;
    }

    @Override
    public boolean equals(Object obj) {
      boolean same = this == obj;
      if (!same && null != obj && ContigStats.class == obj.getClass()) {
        ContigStats rhs = (ContigStats) obj;
        return Objects.equals(allKnownFalsePositives(), rhs.allKnownFalsePositives())
            && Objects.equals(concordance(), rhs.concordance())
            && Objects.equals(contig(), rhs.contig())
            && Objects.equals(correctKnownFalsePositives(), rhs.correctKnownFalsePositives())
            && Objects.equals(falseNegatives(), rhs.falseNegatives())
            && Objects.equals(falsePositives(), rhs.falsePositives())
            && Objects.equals(incorrectPredictions(), rhs.incorrectPredictions())
            && Objects.equals(knownFalsePositives(), rhs.knownFalsePositives())
            && Objects.equals(predictedVariants(), rhs.predictedVariants())
            && Objects.equals(rescuedVariants(), rhs.rescuedVariants())
            && Objects.equals(truePositives(), rhs.truePositives())
            && Objects.equals(trueVariants(), rhs.trueVariants());
      }
      return same;
    }

    public NavigableMap<Integer, Variant> falseNegatives() {
      return falseNegatives;
    }

    public NavigableMap<Integer, Variant> falsePositives() {
      return falsePositives;
    }

    @Override
    public int hashCode() {
      return Objects.hash(
          allKnownFalsePositives(),
          concordance(),
          contig(),
          correctKnownFalsePositives(),
          falseNegatives(),
          falsePositives(),
          incorrectPredictions(),
          knownFalsePositives(),
          predictedVariants(),
          rescuedVariants(),
          truePositives(),
          trueVariants());
    }

    public Multimap<VariantType, Integer> incorrectPredictions() {
      return incorrectPredictions;
    }

    public Optional<NavigableMap<Integer, Variant>> knownFalsePositives() {
      return knownFalsePositives;
    }

    public NavigableMap<Integer, Variant> predictedVariants() {
      return predictedVariants;
    }

    ContigStats rescue(FastaReader.FastaFile reference, int rescueWindowSize) {
      for (Map.Entry<Integer, Variant> entry :
          Maps.newTreeMap(this.falseNegatives).entrySet()) {
        if (this.falseNegatives.containsKey(entry.getKey())) {
          Variant variant = entry.getValue();
          if (!getType.apply(variant).isStructuralVariant()) {
            Optional<SequenceRescuer.RescuedVariants> optional = SequenceRescuer.builder()
                .setContig(contig)
                .setReference(reference)
                .setRescueWindowSize(rescueWindowSize)
                .setTruePositives(truePositives)
                .setFalsePositives(falsePositives)
                .setFalseNegatives(falseNegatives)
                .setGetTypeFunction(getType)
                .build()
                .tryRescue(variant.position());
            if (optional.isPresent()) {
              SequenceRescuer.RescuedVariants rescuedVariants = optional.get();
              NavigableMap<Integer, Variant>
                  newTruePositives = rescuedVariants.newTruePositives(),
                  removeFalsePositives = rescuedVariants.removeFalsePositives();
              for (Integer location : newTruePositives.keySet()) {
                this.falseNegatives.remove(location);
              }
              for (Integer location : removeFalsePositives.keySet()) {
                this.falsePositives.remove(location);
              }
              if (this.rescuedVariants.isPresent()) {
                this.rescuedVariants.get().putAll(newTruePositives);
              } else {
                this.rescuedVariants = Optional.of(newTruePositives);
              }
            }
          }
        }
      }
      return this;
    }

    public Optional<NavigableMap<Integer, Variant>> rescuedVariants() {
      return rescuedVariants;
    }

    ContigStats setKnownFalsePositives(
        NavigableMap<Integer, Variant> knownFalsePositives,
        Iterable<Integer> correctKnownFalsePositiveLocations,
        Iterable<Integer> allKnownFalsePositiveLocations) {
      this.knownFalsePositives =
          Optional.of(knownFalsePositives);
      this.correctKnownFalsePositives =
          Optional.of(filter(knownFalsePositives, correctKnownFalsePositiveLocations));
      this.allKnownFalsePositives =
          Optional.of(filter(knownFalsePositives, allKnownFalsePositiveLocations));
      return this;
    }

    ContigStats setRescuedVariants(NavigableMap<Integer, Variant> rescuedVariants) {
      this.rescuedVariants = Optional.of(rescuedVariants);
      return this;
    }

    @Override
    public String toString() {
      return String.format(
          "allKnownFalsePositives = %s, "
              + "correctKnownFalsePositives = %s, "
              + "knownFalsePositives = %s, "
              + "rescuedVariants = %s, "
              + "concordance = %s, "
              + "contig = %s, "
              + "falseNegatives = %s, "
              + "falsePositives = %s, "
              + "predictedVariants = %s, "
              + "truePositives = %s, "
              + "trueVariants = %s, "
              + "incorrectPredictions = %s",
          toString(allKnownFalsePositives()),
          toString(correctKnownFalsePositives()),
          toString(knownFalsePositives()),
          toString(rescuedVariants()),
          concordance(),
          contig(),
          toString(falseNegatives()),
          toString(falsePositives()),
          toString(predictedVariants()),
          toString(truePositives()),
          toString(trueVariants()),
          incorrectPredictions());
    }

    public NavigableMap<Integer, Variant> truePositives() {
      return truePositives;
    }

    public NavigableMap<Integer, Variant> trueVariants() {
      return trueVariants;
    }
  }

  public static class GenotypeConcordance {

    static GenotypeConcordance create() {
      return new GenotypeConcordance();
    }

    private final Map<VariantType, Map<Genotype, Multiset<Genotype>>> concordance = new HashMap<>();

    private GenotypeConcordance() {
      for (VariantType variantType : VariantType.values()) {
        Map<Genotype, Multiset<Genotype>> map = new HashMap<>();
        for (Genotype genotype : Genotype.values()) {
          map.put(genotype, HashMultiset.<Genotype>create());
        }
        concordance.put(variantType, map);
      }
    }

    @Override
    public boolean equals(Object obj) {
      return this == obj
          || null != obj
          && GenotypeConcordance.class == obj.getClass()
          && Objects.equals(concordance, ((GenotypeConcordance) obj).concordance);
    }

    public int get(VariantType variantType, Genotype trueGenotype, Genotype predictedGenotype) {
      return concordance.get(variantType).get(trueGenotype).count(predictedGenotype);
    }

    @Override
    public int hashCode() {
      return Objects.hashCode(concordance);
    }

    GenotypeConcordance increment(
        VariantType variantType, Genotype trueGenotype, Genotype predictedGenotype) {
      concordance.get(variantType).get(trueGenotype).add(predictedGenotype);
      return this;
    }

    @Override
    public String toString() {
      List<String> parts = new ArrayList<>();
      for (VariantType variantType : VariantType.values()) {
        for (Genotype trueGenotype : Genotype.values()) {
          for (Genotype predictedGenotype : Genotype.values()) {
            int count = get(variantType, trueGenotype, predictedGenotype);
            if (0 < count) {
              parts.add(String.format(
                  "(%s, %s, %s, %d)", variantType, trueGenotype, predictedGenotype, count));
            }
          }
        }
      }
      return Joiner.on(", ").join(parts);
    }
  }

  private static final Function<Variant, Integer>
      GET_POSITION =
          new Function<Variant, Integer>() {
            @Override public Integer apply(Variant variant) {
              return variant.position();
            }
          };

  private static final Function<Iterable<Variant>, Multimap<String, Variant>>
      INDEX_BY_CONTIG =
          new Function<Iterable<Variant>, Multimap<String, Variant>>() {
            @Override public Multimap<String, Variant> apply(Iterable<Variant> variants) {
              return Multimaps.index(variants, Variant.CONTIG);
            }
          };

  private static final Function<Iterable<Variant>, NavigableMap<Integer, Variant>>
      INDEX_BY_POSITION = uniqueIndex(GET_POSITION);

  private static final Function<Multimap<String, ?>, Set<String>>
      KEY_SET = keySetFunction();

  public static Builder builder() {
    return new Builder();
  }

  private static <X extends Comparable<? super X>> NavigableSet<X>
      difference(Set<? extends X> lhs, Set<? extends X> rhs) {
    NavigableSet<X> difference = new TreeSet<>();
    difference.addAll(lhs);
    difference.removeAll(rhs);
    return difference;
  }

  @SafeVarargs
  private static <X extends Comparable<? super X>> NavigableSet<X>
      intersection(Set<? extends X>... sets) {
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

  private static <X> Function<Multimap<X, ?>, Set<X>> keySetFunction() {
    return new Function<Multimap<X, ?>, Set<X>>() {
          @Override public Set<X> apply(Multimap<X, ?> multimap) {
            return multimap.keySet();
          }
        };
  }

  @SafeVarargs
  private static <X extends Comparable<? super X>> NavigableSet<X>
      union(Set<? extends X>... sets) {
    NavigableSet<X> union = new TreeSet<>();
    for (Set<? extends X> set : sets) {
      union.addAll(set);
    }
    return union;
  }

  private static <X extends Comparable<? super X>, Y> Function<Iterable<Y>, NavigableMap<X, Y>>
      uniqueIndex(final Function<? super Y, X> function) {
    return new Function<Iterable<Y>, NavigableMap<X, Y>>() {
          @Override public NavigableMap<X, Y> apply(Iterable<Y> iterable) {
            NavigableMap<X, Y> result = new TreeMap<>();
            for (Y object : iterable) {
              result.put(function.apply(object), object);
            }
            return result;
          }
        };
  }

  private final Function<Variant, VariantType> getType;
  private final int maxSvBreakpointDistance;
  private final int maxVariantLengthDifference;
  private final Optional<FastaReader.FastaFile> reference;
  private final int rescueWindowSize;

  private VariantEvaluator(
      Optional<FastaReader.FastaFile> reference,
      int maxSvBreakpointDistance,
      int maxVariantLengthDifference,
      int rescueWindowSize,
      Function<Variant, VariantType> getType) {
    this.reference = reference;
    this.maxSvBreakpointDistance = maxSvBreakpointDistance;
    this.maxVariantLengthDifference = maxVariantLengthDifference;
    this.rescueWindowSize = rescueWindowSize;
    this.getType = getType;
  }

  public Map<String, ContigStats> evaluate(
      Iterable<Variant> trueVariants,
      Iterable<Variant> predictedVariants) {
    return evaluate(trueVariants, predictedVariants, Optional.<Iterable<Variant>>absent());
  }

  public Map<String, ContigStats> evaluate(
      Iterable<Variant> trueVariants,
      Iterable<Variant> predictedVariants,
      Iterable<Variant> knownFalsePositives) {
    return evaluate(trueVariants, predictedVariants, Optional.of(knownFalsePositives));
  }

  Map<String, ContigStats> evaluate(
      Iterable<Variant> trueVariants,
      Iterable<Variant> predictedVariants,
      Optional<Iterable<Variant>> knownFalsePositives) {
    Multimap<String, Variant>
        trueVars = INDEX_BY_CONTIG.apply(trueVariants),
        predVars = INDEX_BY_CONTIG.apply(predictedVariants);
    Optional<Multimap<String, Variant>>
        knownFps = knownFalsePositives.transform(INDEX_BY_CONTIG);
    ImmutableMap.Builder<String, ContigStats> result = ImmutableMap.builder();
    for (final String contig : union(
        KEY_SET.apply(trueVars),
        KEY_SET.apply(predVars),
        knownFps.transform(KEY_SET).or(Collections.<String>emptySet()))) {
      Function<Multimap<String, Variant>, NavigableMap<Integer, Variant>>
          indexByPosition = Functions.compose(
              INDEX_BY_POSITION,
              new Function<Multimap<String, Variant>, Iterable<Variant>>() {
                @Override public Iterable<Variant> apply(
                    Multimap<String, Variant> multimap) {
                  return multimap.get(contig);
                }
              });
      result.put(
          contig,
          evaluate(
              contig,
              indexByPosition.apply(trueVars),
              indexByPosition.apply(predVars),
              knownFps.transform(indexByPosition)));
    }
    return result.build();
  }

  private ContigStats evaluate(
      String contig,
      NavigableMap<Integer, Variant> trueVariants,
      NavigableMap<Integer, Variant> predictedVariants,
      Optional<NavigableMap<Integer, Variant>> knownFalsePositives) {
    Set<Integer>
        trueVariantLocations = trueVariants.keySet(),
        predictedVariantLocations = predictedVariants.keySet(),
        truePositiveLocations = new HashSet<>(),
        falsePositiveLocations = difference(predictedVariantLocations, trueVariantLocations),
        falseNegativeLocations = difference(trueVariantLocations, predictedVariantLocations),
        incorrectPredictions = new HashSet<>(),
        correctKnownFalsePositiveLocations = new HashSet<>(),
        allKnownFalsePositiveLocations = new HashSet<>();
    GenotypeConcordance concordance = GenotypeConcordance.create();
    boolean knownFalsePositivesPresent = knownFalsePositives.isPresent();
    for (Integer location : intersection(predictedVariantLocations, trueVariantLocations)) {
      Variant
          trueVariant = trueVariants.get(location),
          predictedVariant = predictedVariants.get(location);
      VariantType trueVariantType = getType.apply(trueVariant);
      if (trueVariantType == getType.apply(predictedVariant)
          && (trueVariantType.isStructuralVariant() || trueVariant.alternateBases()
              .equals(predictedVariant.alternateBases()))) {
        truePositiveLocations.add(location);
        concordance.increment(
            trueVariantType,
            getGenotype(trueVariant),
            getGenotype(predictedVariant));
      } else {
        incorrectPredictions.add(location);
      }
    }
    if (knownFalsePositivesPresent) {
      NavigableMap<Integer, Variant> knownFp = knownFalsePositives.get();
      for (Integer location :
          intersection(predictedVariantLocations, knownFp.keySet())) {
        Variant knownFalsePositive = knownFp.get(location);
        if (Objects.equals(
            knownFalsePositive.referenceBases(),
            predictedVariants.get(location).referenceBases())) {
          correctKnownFalsePositiveLocations.add(location);
        }
        allKnownFalsePositiveLocations.add(location);
      }
    }
    for (Integer location : difference(trueVariantLocations, predictedVariantLocations)) {
      Variant trueVariant = trueVariants.get(location);
      VariantType trueVariantType = getType.apply(trueVariant);
      if (trueVariantType.isStructuralVariant()) {
        Optional<Integer> structuralMatch =
            structuralMatch(trueVariant, predictedVariants).transform(GET_POSITION);
        if (structuralMatch.isPresent()) {
          Integer match = structuralMatch.get();
          if (falsePositiveLocations.contains(match)) {
            truePositiveLocations.add(location);
            falsePositiveLocations.remove(match);
            falseNegativeLocations.remove(location);
            concordance.increment(
                trueVariantType,
                getGenotype(trueVariant),
                getGenotype(predictedVariants.get(match)));
          }
        }
      }
    }
    falsePositiveLocations.addAll(incorrectPredictions);
    falseNegativeLocations.addAll(incorrectPredictions);
    ContigStats variantStats = ContigStats.create(contig, trueVariants, predictedVariants,
        truePositiveLocations, falsePositiveLocations, falseNegativeLocations, incorrectPredictions,
        concordance, getType);
    if (reference.isPresent()) {
      variantStats.rescue(reference.get(), rescueWindowSize);
    }
    if (knownFalsePositivesPresent) {
      variantStats.setKnownFalsePositives(
          knownFalsePositives.get(),
          correctKnownFalsePositiveLocations,
          allKnownFalsePositiveLocations);
    }
    return variantStats;
  }

  private Optional<Variant> structuralMatch(
      final Variant trueVariant,
      NavigableMap<Integer, Variant> predictedVariants) {
    final int trueVariantPosition = trueVariant.position();
    Collection<Variant> candidates = FluentIterable
        .from(predictedVariants
            .subMap(
                trueVariantPosition - maxSvBreakpointDistance, true,
                trueVariantPosition + maxSvBreakpointDistance, true)
            .values())
        .filter(
            new Predicate<Variant>() {

              private final VariantType trueVariantType = getType.apply(trueVariant);

              @Override public boolean apply(Variant predictedVariant) {
                if (trueVariantType == getType.apply(predictedVariant)) {
                  int trueReferenceLength = trueVariant.referenceBases().length(),
                      predictedReferenceLength = predictedVariant.referenceBases().length();
                  for (String trueAlt : trueVariant.alternateBases()) {
                    int trueGain = trueReferenceLength - trueAlt.length();
                    for (String predictedAlt : predictedVariant.alternateBases()) {
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
    return candidates.isEmpty()
        ? Optional.<Variant>absent()
        : Optional.of(Ordering
            .from(
                new Comparator<Variant>() {
                  @Override public int compare(Variant lhs, Variant rhs) {
                    return Long.compare(
                        Math.abs(trueVariantPosition - lhs.position()),
                        Math.abs(trueVariantPosition - rhs.position()));
                  }
                })
            .min(candidates));
  }

  private static Genotype getGenotype(Variant variant) {
    return Iterables.getOnlyElement(Genotype.getGenotypes(variant));
  }
}

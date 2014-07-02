package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.AbstractSequentialIterator;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;

import edu.berkeley.cs.amplab.smash4j.fasta.FastaReader;
import edu.berkeley.cs.amplab.smash4j.fasta.FastaReader.FastaFile;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Objects;
import java.util.TreeMap;

public class SequenceRescuer {

  public static class Builder {

    private String contig;
    private NavigableMap<Integer, Variant> falseNegatives;
    private NavigableMap<Integer, Variant> falsePositives;
    private Function<Variant, VariantType> getType;
    private FastaReader.FastaFile reference;
    private int rescueWindowSize;
    private NavigableMap<Integer, Variant> truePositives;

    public SequenceRescuer build() {
      return new SequenceRescuer(
          Preconditions.checkNotNull(contig),
          Preconditions.checkNotNull(truePositives),
          Preconditions.checkNotNull(falsePositives),
          Preconditions.checkNotNull(falseNegatives),
          Preconditions.checkNotNull(reference),
          Preconditions.checkNotNull(Window.Factory.builder()
              .setTruePositives(truePositives)
              .setFalsePositives(falsePositives)
              .setFalseNegatives(falseNegatives)
              .setSize(rescueWindowSize)
              .build()),
              getType);
    }

    public Builder setContig(String contig) {
      this.contig = contig;
      return this;
    }

    public Builder setFalseNegatives(NavigableMap<Integer, Variant> falseNegatives) {
      this.falseNegatives = falseNegatives;
      return this;
    }

    public Builder setFalsePositives(NavigableMap<Integer, Variant> falsePositives) {
      this.falsePositives = falsePositives;
      return this;
    }

    public Builder setGetTypeFunction(Function<Variant, VariantType> getType) {
      this.getType = getType;
      return this;
    }

    public Builder setReference(FastaReader.FastaFile reference) {
      this.reference = reference;
      return this;
    }

    public Builder setRescueWindowSize(int rescueWindowSize) {
      this.rescueWindowSize = rescueWindowSize;
      return this;
    }

    public Builder setTruePositives(NavigableMap<Integer, Variant> truePositives) {
      this.truePositives = truePositives;
      return this;
    }
  }

  public static class RescuedVariants {

    static RescuedVariants create(
        NavigableMap<Integer, Variant> newTruePositives,
        NavigableMap<Integer, Variant> removeFalsePositives) {
      return new RescuedVariants(newTruePositives, removeFalsePositives);
    }

    private final NavigableMap<Integer, Variant> newTruePositives, removeFalsePositives;

    private RescuedVariants(
        NavigableMap<Integer, Variant> newTruePositives,
        NavigableMap<Integer, Variant> removeFalsePositives) {
      this.newTruePositives = newTruePositives;
      this.removeFalsePositives = removeFalsePositives;
    }

    @Override
    public boolean equals(Object obj) {
      boolean same = this == obj;
      if (!same && null != obj && RescuedVariants.class == obj.getClass()) {
        RescuedVariants rhs = (RescuedVariants) obj;
        return Objects.equals(newTruePositives(), rhs.newTruePositives())
            && Objects.equals(removeFalsePositives(), rhs.removeFalsePositives());
      }
      return same;
    }

    @Override
    public int hashCode() {
      return Objects.hash(newTruePositives(), removeFalsePositives());
    }

    public NavigableMap<Integer, Variant> newTruePositives() {
      return newTruePositives;
    }

    public NavigableMap<Integer, Variant> removeFalsePositives() {
      return removeFalsePositives;
    }
  }

  static class Window {

    static class Factory {

      static class Builder {

        private static <X> Function<X, X> fix(final Function<? super X, ? extends X> function) {
          return new Function<X, X>() {
                @Override public X apply(X initial) {
                  return Iterators.getLast(
                      new AbstractSequentialIterator<X>(initial) {
                        @Override protected X computeNext(X previous) {
                          X next = function.apply(previous);
                          return Objects.equals(previous, next) ? null : next;
                        }
                      });
                }
              };
        }

        static Function<Window, Window> windowEnlarger(
            final NavigableMap<Integer, Variant> variants) {
          return new Function<Window, Window>() {
                @Override public Window apply(Window window) {
                  final int
                      lowerBound = window.lowerBound(),
                      upperBound = window.upperBound();
                  return Window.create(
                      getChoppedVariant(variants, lowerBound, LOWEST_START)
                          .transform(GET_START)
                          .or(lowerBound),
                      getChoppedVariant(variants, upperBound - 1, HIGHEST_END)
                          .transform(
                              new Function<Variant, Integer>() {
                                @Override public Integer apply(Variant variant) {
                                  int end = GET_END.apply(variant);
                                  return upperBound == end ? end + 1 : end;
                                }
                              })
                          .or(upperBound));
                }
              };
        }

        private NavigableMap<Integer, Variant> falseNegatives, falsePositives, truePositives;
        private int size;

        Factory build() {
          return new Factory(fix(Functions.compose(Functions.compose(windowEnlarger(truePositives),
              windowEnlarger(falsePositives)), windowEnlarger(falseNegatives))), size);
        }

        Builder setFalseNegatives(NavigableMap<Integer, Variant> falseNegatives) {
          this.falseNegatives = falseNegatives;
          return this;
        }

        Builder setFalsePositives(NavigableMap<Integer, Variant> falsePositives) {
          this.falsePositives = falsePositives;
          return this;
        }

        Builder setSize(int size) {
          this.size = size;
          return this;
        }

        Builder setTruePositives(NavigableMap<Integer, Variant> truePositives) {
          this.truePositives = truePositives;
          return this;
        }
      }

      private static final int
          WINDOW_SIZE_LIMIT = 5000;

      static Builder builder() {
        return new Builder();
      }

      private final int size;
      private final Function<Window, Window> windowEnlarger;

      private Factory(Function<Window, Window> windowEnlarger, int size) {
        this.windowEnlarger = windowEnlarger;
        this.size = size;
      }

      Optional<Window> createWindow(int location) {
        Window window = windowEnlarger.apply(Window.create(location - size, location + size));
        return WINDOW_SIZE_LIMIT < window.upperBound() - window.lowerBound()
            ? Optional.<Window>absent()
            : Optional.of(window);
      }
    }

    static Window create(int lowerBound, int upperBound) {
      return new Window(lowerBound, upperBound);
    }

    private final int lowerBound, upperBound;

    private Window(int lowerBound, int upperBound) {
      this.lowerBound = lowerBound;
      this.upperBound = upperBound;
    }

    @Override public boolean equals(Object obj) {
      boolean same = this == obj;
      if (!same && null != obj && Window.class == obj.getClass()) {
        Window rhs = (Window) obj;
        return lowerBound() == rhs.lowerBound() && upperBound() == rhs.upperBound();
      }
      return same;
    }

    @Override public int hashCode() {
      return Objects.hash(lowerBound(), upperBound());
    }

    int lowerBound() {
      return lowerBound;
    }

    <X> NavigableMap<Integer, X> restrict(NavigableMap<Integer, X> map) {
      return map.subMap(lowerBound(), true, upperBound(), false);
    }

    @Override public String toString() {
      return String.format("[%d, %d]", lowerBound(), upperBound());
    }

    int upperBound() {
      return upperBound;
    }
  }

  private static final Function<Iterable<? extends Variant>, List<List<Variant>>>
      GET_OVERLAPS = grouper(
          new Function<Variant, Predicate<Variant>>() {
            @Override public Predicate<Variant> apply(final Variant candidate) {
              return strictlyOverlaps(candidate.position());
            }
          });

  private static final Function<Variant, Integer>
      GET_START =
          new Function<Variant, Integer>() {
            @Override public Integer apply(Variant variant) {
              return variant.position();
            }
          },
      GET_END =
          new Function<Variant, Integer>() {
            @Override public Integer apply(Variant variant) {
              return GET_START.apply(variant) + variant.referenceBases().length();
            }
          };

  static final Ordering<Variant>
      LOWEST_START = Ordering.natural().onResultOf(GET_START),
      HIGHEST_END = Ordering.natural().onResultOf(GET_END).reverse();

  private static final int
      WINDOW_MAX_OVERLAPPING = 16,
      WINDOW_VARIANT_LOOKBACK_SIZE = 50;

  private static StringBuilder addRefBasesUntil(StringBuilder chunks,
      FastaReader.FastaFile reference, String contig, int begin, int end) {
    return chunks.append(getRefBases(reference, contig, begin, end));
  }

  private static Optional<List<Variant>> addTruePosToQueue(
      List<Variant> queue, List<Variant> truePositives) {
    List<Variant> newQueue = new ArrayList<>();
    newQueue.addAll(queue);
    for (Variant truePositive : truePositives) {
      for (Variant variant : queue) {
        if (strictlyOverlaps(truePositive.position()).apply(variant)) {
          return Optional.absent();
        }
      }
      newQueue.add(truePositive);
    }
    Collections.sort(newQueue, Ordering.natural().onResultOf(GET_START));
    return Optional.of(newQueue);
  }

  public static Builder builder() {
    return new Builder();
  }

  private static <X extends Comparable<? super X>, Y> NavigableMap<X, Y>
      filter(Map<? extends X, ? extends Y> map, Predicate<? super Y> predicate) {
    NavigableMap<X, Y> filtered = new TreeMap<>();
    for (Map.Entry<? extends X, ? extends Y> entry : map.entrySet()) {
      Y value = entry.getValue();
      if (predicate.apply(value)) {
        filtered.put(entry.getKey(), value);
      }
    }
    return filtered;
  }

  static Optional<Variant> getChoppedVariant(
      NavigableMap<Integer, Variant> variants,
      int location, Ordering<Variant> ordering) {
    Collection<Variant> collection = FluentIterable
        .from(variants
            .subMap(location - WINDOW_VARIANT_LOOKBACK_SIZE, true, location, true)
            .values())
        .filter(overlaps(location))
        .toList();
    return collection.isEmpty()
        ? Optional.<Variant>absent()
        : Optional.of(ordering.min(collection));
  }

  private static String getRefBases(
      FastaReader.FastaFile reference, String contig, int begin, int end) {
    return reference.get(contig, begin - 1, end - 1);
  }

  private static List<List<Variant>> getRestOfPath(
      final List<Variant> chosenSoFar,
      final List<List<Variant>> remainingChoices) {
    return remainingChoices.isEmpty()
        ? Collections.singletonList(chosenSoFar)
        : ImmutableList.copyOf(Iterables.concat(FluentIterable
            .from(remainingChoices.get(0))
            .filter(
                new Predicate<Variant>() {
                  @Override public boolean apply(Variant choice) {
                    return !Iterables.any(chosenSoFar, strictlyOverlaps(choice.position()));
                  }
                })
            .transform(
                new Function<Variant, List<List<Variant>>>() {
                  @Override public List<List<Variant>> apply(Variant choice) {
                    return getRestOfPath(
                        ImmutableList.<Variant>builder()
                            .addAll(chosenSoFar)
                            .add(choice)
                            .build(),
                        remainingChoices.subList(1, remainingChoices.size()));
                  }
                })));
  }

  static Optional<String> getSequence(
      FastaReader.FastaFile reference, String contig, Window window, List<Variant> variants) {
    StringBuilder builder = new StringBuilder();
    int homOffset = window.lowerBound();
    for (Variant variant : variants) {
      String referenceBases = variant.referenceBases();
      int position = variant.position(), next = referenceBases.length() + position;
      if (!Objects.equals(referenceBases, getRefBases(reference, contig, position, next))) {
        return Optional.absent();
      }
      addRefBasesUntil(builder, reference, contig, homOffset, position).append(
          variant.alternateBases().get(0));
      homOffset = next;
    }
    return Optional.of(
        addRefBasesUntil(builder, reference, contig, homOffset, window.upperBound()).toString());
  }

  private static <X> Function<Iterable<? extends X>, List<List<X>>> grouper(
      final Function<? super X, ? extends Predicate<? super X>> function) {
    return new Function<Iterable<? extends X>, List<List<X>>>() {
          @Override public List<List<X>> apply(final Iterable<? extends X> iterable) {
            return ImmutableList.copyOf(
                new AbstractIterator<List<X>>() {

                  private final PeekingIterator<X>
                      iterator = Iterators.peekingIterator(iterable.iterator());

                  @Override protected List<X> computeNext() {
                    if (iterator.hasNext()) {
                      List<X> list = new ArrayList<>();
                      for (
                          list.add(iterator.next());
                          iterator.hasNext()
                              && Iterables.any(list, function.apply(iterator.peek()));
                          list.add(iterator.next()));
                      return list;
                    }
                    return endOfData();
                  }
                });
          }
        };
  }

  private final Predicate<Variant>
      isNonSnp = Predicates.not(
          new Predicate<Variant>() {
            @Override public boolean apply(Variant variant) {
              return getType.apply(variant).isSnp();
            }
          }),
      isNotStructuralVariant = Predicates.not(
          new Predicate<Variant>() {
            @Override public boolean apply(Variant variant) {
              return getType.apply(variant).isStructuralVariant();
            }
          });

  private static List<Integer> losses(Variant variant) {
    int referenceBasesLength = variant.referenceBases().length();
    ImmutableList.Builder<Integer> losses = ImmutableList.builder();
    for (String alternateBases : variant.alternateBases()) {
      losses.add(Math.max(0, alternateBases.length() - referenceBasesLength));
    }
    return losses.build();
  }

  private static Predicate<Variant> overlaps(final int location) {
    return new Predicate<Variant>() {
          @Override public boolean apply(Variant variant) {
            return GET_START.apply(variant) <= location
                && location < GET_END.apply(variant);
          }
        };
  }

  private static Predicate<Variant> overlapsAllele(final int location) {
    return new Predicate<Variant>() {
          @Override public boolean apply(Variant variant) {
            int position = variant.position();
            for (Integer loss : losses(variant)) {
              if (location - loss <= position && position <= location) {
                return true;
              }
            }
            return false;
          }
        };
  }

  private static Predicate<Variant> strictlyOverlaps(final int location) {
    return new Predicate<Variant>() {
          @Override public boolean apply(Variant variant) {
            int position = variant.position();
            return position <= location
                && location < position + variant.referenceBases().length();
          }
        };
  }

  private static <X extends Comparable<? super X>, Y> NavigableMap<X, Y> uniqueIndex(
      Iterable<? extends Y> iterable, Function<? super Y, ? extends X> function) {
    NavigableMap<X, Y> index = new TreeMap<>();
    for (Y object : iterable) {
      index.put(function.apply(object), object);
    }
    return index;
  }

  private final String contig;
  private final NavigableMap<Integer, Variant> falseNegatives;
  private final NavigableMap<Integer, Variant> falsePositives;
  private final Function<Variant, VariantType> getType;
  private final FastaReader.FastaFile reference;
  private final NavigableMap<Integer, Variant> truePositives;
  private final Window.Factory windowFactory;

  private SequenceRescuer(
      String contig,
      NavigableMap<Integer, Variant> truePositives,
      NavigableMap<Integer, Variant> falsePositives,
      NavigableMap<Integer, Variant> falseNegatives,
      FastaFile reference,
      Window.Factory windowFactory,
      Function<Variant, VariantType> getType) {
    this.contig = contig;
    this.truePositives = truePositives;
    this.falsePositives = falsePositives;
    this.falseNegatives = falseNegatives;
    this.reference = reference;
    this.windowFactory = windowFactory;
    this.getType = getType;
  }

  private List<Variant> extractRangeAndFilter(
      NavigableMap<Integer, Variant> variants, Window window, int location) {
    NavigableMap<Integer, Variant>
        result = filter(window.restrict(variants), isNotStructuralVariant);
    Variant variant = variants.get(location);
    if (null == variant) {
      return ImmutableList.copyOf(result.values());
    }
    Predicate<Variant> predicate = Predicates.or(
        Predicates.compose(Predicates.equalTo(location), GET_START),
        Predicates.not(Predicates.or(
            overlapsAllele(location),
            overlapsAllele(location + Ordering.natural().max(losses(variant))))));
    return ImmutableList.copyOf(filter(result, predicate).values());
  }

  private List<List<Variant>> extractVariantQueues(
      NavigableMap<Integer, Variant> variants, Window window, int location) {
    List<Variant>
        variantsInWindow = extractRangeAndFilter(variants, window, location);
    return variantsInWindow.isEmpty()
        ? Collections.<List<Variant>>emptyList()
        : getRestOfPath(new ArrayList<Variant>(), GET_OVERLAPS.apply(variantsInWindow));
  }

  public Optional<RescuedVariants> tryRescue(final int position) {
    return windowFactory.createWindow(position)
        .transform(
            new Function<Window, Optional<RescuedVariants>>() {
              @Override public Optional<RescuedVariants> apply(Window window) {
                return tryRescue(position, window);
              }
            })
        .or(Optional.<RescuedVariants>absent());
  }

  private Optional<RescuedVariants> tryRescue(int location, final Window window) {
    List<List<Variant>>
        falseNegativesQueue = extractVariantQueues(this.falseNegatives, window, location),
        falsePositivesQueue = extractVariantQueues(this.falsePositives, window, location);
    final List<Variant> truePositives =
        extractRangeAndFilter(this.truePositives, window, location);
    if (falseNegativesQueue.isEmpty()
        || falsePositivesQueue.isEmpty()
        || WINDOW_MAX_OVERLAPPING < falseNegativesQueue.size() * falsePositivesQueue.size()) {
      return Optional.absent();
    }
    for (final List<Variant> newTruePositives : falseNegativesQueue) {
      for (final List<Variant> removeFalsePositives : falsePositivesQueue) {
        if (Iterables.any(Iterables.concat(newTruePositives, removeFalsePositives), isNonSnp)) {
          return addTruePosToQueue(newTruePositives, truePositives)
              .transform(
                  new Function<List<Variant>, Optional<RescuedVariants>>() {
                    @Override public Optional<RescuedVariants> apply(
                        final List<Variant> falseNegatives) {
                      return addTruePosToQueue(removeFalsePositives, truePositives)
                          .transform(
                              new Function<List<Variant>, Optional<RescuedVariants>>() {
                                @Override public Optional<RescuedVariants> apply(
                                    List<Variant> falsePositives) {
                                  Optional<String>
                                      falseNegativesSequence =
                                          getSequence(reference, contig, window, falseNegatives),
                                      falsePositivesSequence =
                                          getSequence(reference, contig, window, falsePositives);
                                  return falseNegativesSequence.isPresent()
                                      && falsePositivesSequence.isPresent()
                                      && Objects.equals(
                                              falseNegativesSequence.get(),
                                              falsePositivesSequence.get())
                                          ? Optional.of(RescuedVariants.create(
                                              uniqueIndex(newTruePositives, GET_START),
                                              uniqueIndex(removeFalsePositives, GET_START)))
                                          : Optional.<RescuedVariants>absent();
                                }
                              })
                          .or(Optional.<RescuedVariants>absent());
                    }
                  })
              .or(Optional.<RescuedVariants>absent());
        }
      }
    }
    return Optional.absent();
  }
}

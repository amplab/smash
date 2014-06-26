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

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.fastaparser.FastaReader.Callback.FastaFile;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.smash4j.VariantEvaluator.VariantType;

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
    private NavigableMap<Integer, VariantProto> falseNegatives;
    private NavigableMap<Integer, VariantProto> falsePositives;
    private FastaReader.Callback.FastaFile reference;
    private int rescueWindowSize;
    private NavigableMap<Integer, VariantProto> truePositives;

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
              .build()));
    }

    public Builder setContig(String contig) {
      this.contig = contig;
      return this;
    }

    public Builder setFalseNegatives(NavigableMap<Integer, VariantProto> falseNegatives) {
      this.falseNegatives = falseNegatives;
      return this;
    }

    public Builder setFalsePositives(NavigableMap<Integer, VariantProto> falsePositives) {
      this.falsePositives = falsePositives;
      return this;
    }

    public Builder setReference(FastaReader.Callback.FastaFile reference) {
      this.reference = reference;
      return this;
    }

    public Builder setRescueWindowSize(int rescueWindowSize) {
      this.rescueWindowSize = rescueWindowSize;
      return this;
    }

    public Builder setTruePositives(NavigableMap<Integer, VariantProto> truePositives) {
      this.truePositives = truePositives;
      return this;
    }
  }

  public static class RescuedVariants {

    static RescuedVariants create(
        NavigableMap<Integer, VariantProto> newTruePositives,
        NavigableMap<Integer, VariantProto> removeFalsePositives) {
      return new RescuedVariants(newTruePositives, removeFalsePositives);
    }

    private final NavigableMap<Integer, VariantProto> newTruePositives, removeFalsePositives;

    private RescuedVariants(
        NavigableMap<Integer, VariantProto> newTruePositives,
        NavigableMap<Integer, VariantProto> removeFalsePositives) {
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

    public NavigableMap<Integer, VariantProto> newTruePositives() {
      return newTruePositives;
    }

    public NavigableMap<Integer, VariantProto> removeFalsePositives() {
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
            final NavigableMap<Integer, VariantProto> variants) {
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
                              new Function<VariantProto, Integer>() {
                                @Override public Integer apply(VariantProto variant) {
                                  int end = GET_END.apply(variant);
                                  return upperBound == end ? end + 1 : end;
                                }
                              })
                          .or(upperBound));
                }
              };
        }

        private NavigableMap<Integer, VariantProto> falseNegatives, falsePositives, truePositives;
        private int size;

        Factory build() {
          return new Factory(fix(Functions.compose(Functions.compose(windowEnlarger(truePositives),
              windowEnlarger(falsePositives)), windowEnlarger(falseNegatives))), size);
        }

        Builder setFalseNegatives(NavigableMap<Integer, VariantProto> falseNegatives) {
          this.falseNegatives = falseNegatives;
          return this;
        }

        Builder setFalsePositives(NavigableMap<Integer, VariantProto> falsePositives) {
          this.falsePositives = falsePositives;
          return this;
        }

        Builder setSize(int size) {
          this.size = size;
          return this;
        }

        Builder setTruePositives(NavigableMap<Integer, VariantProto> truePositives) {
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

  private static final Function<Iterable<? extends VariantProto>, List<List<VariantProto>>>
      GET_OVERLAPS = grouper(
          new Function<VariantProto, Predicate<VariantProto>>() {
            @Override public Predicate<VariantProto> apply(final VariantProto candidate) {
              return strictlyOverlaps(candidate.getPosition());
            }
          });

  private static final Function<VariantProto, Integer>
      GET_START =
          new Function<VariantProto, Integer>() {
            @Override public Integer apply(VariantProto variant) {
              return variant.getPosition();
            }
          },
      GET_END =
          new Function<VariantProto, Integer>() {
            @Override public Integer apply(VariantProto variant) {
              return GET_START.apply(variant) + variant.getReferenceBases().length();
            }
          };

  private static final Predicate<VariantProto>
      IS_NON_SNP = Predicates.not(Predicates.compose(
          new Predicate<VariantType>() {
            @Override public boolean apply(VariantType type) {
              return type.isSnp();
            }
          },
          VariantEvaluator.VariantType.GET_TYPE)),
      NOT_SV = Predicates.not(Predicates.compose(
          new Predicate<VariantType>() {
            @Override public boolean apply(VariantType type) {
              return type.isStructuralVariant();
            }
          },
          VariantEvaluator.VariantType.GET_TYPE));

  static final Ordering<VariantProto>
      LOWEST_START = Ordering.natural().onResultOf(GET_START),
      HIGHEST_END = Ordering.natural().onResultOf(GET_END).reverse();

  private static final int
      WINDOW_MAX_OVERLAPPING = 16,
      WINDOW_VARIANT_LOOKBACK_SIZE = 50;

  private static Optional<List<VariantProto>> addTruePosToQueue(
      List<VariantProto> queue, List<VariantProto> truePositives) {
    List<VariantProto> newQueue = new ArrayList<>();
    newQueue.addAll(queue);
    for (VariantProto truePositive : truePositives) {
      for (VariantProto variant : queue) {
        if (strictlyOverlaps(truePositive.getPosition()).apply(variant)) {
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

  private static List<VariantProto> extractRangeAndFilter(
      NavigableMap<Integer, VariantProto> variants, Window window, int location) {
    NavigableMap<Integer, VariantProto> result = filter(window.restrict(variants), NOT_SV);
    VariantProto variant = variants.get(location);
    if (null == variant) {
      return ImmutableList.copyOf(result.values());
    }
    Predicate<VariantProto> predicate = Predicates.or(
        Predicates.compose(Predicates.equalTo(location), GET_START),
        Predicates.not(Predicates.or(
            overlapsAllele(location),
            overlapsAllele(location + Ordering.natural().max(losses(variant))))));
    return ImmutableList.copyOf(filter(result, predicate).values());
  }

  private static List<List<VariantProto>> extractVariantQueues(
      NavigableMap<Integer, VariantProto> variants, Window window, int location) {
    List<VariantProto> variantsInWindow = extractRangeAndFilter(variants, window, location);
    return variantsInWindow.isEmpty()
        ? Collections.<List<VariantProto>>emptyList()
        : getRestOfPath(new ArrayList<VariantProto>(), GET_OVERLAPS.apply(variantsInWindow));
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

  static Optional<VariantProto> getChoppedVariant(
      NavigableMap<Integer, VariantProto> variants,
      int location, Ordering<VariantProto> ordering) {
    Collection<VariantProto> collection = FluentIterable
        .from(variants
            .subMap(location - WINDOW_VARIANT_LOOKBACK_SIZE, true, location, true)
            .values())
        .filter(overlaps(location))
        .toList();
    return collection.isEmpty()
        ? Optional.<VariantProto>absent()
        : Optional.of(ordering.min(collection));
  }

  private static List<List<VariantProto>> getRestOfPath(
      final List<VariantProto> chosenSoFar,
      final List<List<VariantProto>> remainingChoices) {
    return remainingChoices.isEmpty()
        ? Collections.singletonList(chosenSoFar)
        : ImmutableList.copyOf(Iterables.concat(FluentIterable
            .from(remainingChoices.get(0))
            .filter(
                new Predicate<VariantProto>() {
                  @Override public boolean apply(VariantProto choice) {
                    return !Iterables.any(chosenSoFar, strictlyOverlaps(choice.getPosition()));
                  }
                })
            .transform(
                new Function<VariantProto, List<List<VariantProto>>>() {
                  @Override public List<List<VariantProto>> apply(VariantProto choice) {
                    return getRestOfPath(
                        ImmutableList.<VariantProto>builder()
                            .addAll(chosenSoFar)
                            .add(choice)
                            .build(),
                        remainingChoices.subList(1, remainingChoices.size()));
                  }
                })));
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

  private static List<Integer> losses(VariantProto variant) {
    int referenceBasesLength = variant.getReferenceBases().length();
    ImmutableList.Builder<Integer> losses = ImmutableList.builder();
    for (String alternateBases : variant.getAlternateBasesList()) {
      losses.add(Math.max(0, alternateBases.length() - referenceBasesLength));
    }
    return losses.build();
  }

  private static Predicate<VariantProto> overlaps(final int location) {
    return new Predicate<VariantProto>() {
          @Override public boolean apply(VariantProto variant) {
            return GET_START.apply(variant) <= location
                && location < GET_END.apply(variant);
          }
        };
  }

  private static Predicate<VariantProto> overlapsAllele(final int location) {
    return new Predicate<VariantProto>() {
          @Override public boolean apply(VariantProto variant) {
            int position = variant.getPosition();
            for (Integer loss : losses(variant)) {
              if (location - loss <= position && position <= location) {
                return true;
              }
            }
            return false;
          }
        };
  }

  private static Predicate<VariantProto> strictlyOverlaps(final int location) {
    return new Predicate<VariantProto>() {
          @Override public boolean apply(VariantProto variant) {
            int position = variant.getPosition();
            return position <= location
                && location < position + variant.getReferenceBases().length();
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
  private final NavigableMap<Integer, VariantProto> falseNegatives;
  private final NavigableMap<Integer, VariantProto> falsePositives;
  private final FastaReader.Callback.FastaFile reference;
  private final NavigableMap<Integer, VariantProto> truePositives;

  private final Window.Factory windowFactory;

  private SequenceRescuer(
      String contig,
      NavigableMap<Integer, VariantProto> truePositives,
      NavigableMap<Integer, VariantProto> falsePositives,
      NavigableMap<Integer, VariantProto> falseNegatives,
      FastaFile reference,
      Window.Factory windowFactory) {
    this.contig = contig;
    this.truePositives = truePositives;
    this.falsePositives = falsePositives;
    this.falseNegatives = falseNegatives;
    this.reference = reference;
    this.windowFactory = windowFactory;
  }

  private static StringBuilder addRefBasesUntil(StringBuilder chunks,
      FastaReader.Callback.FastaFile reference, String contig, int begin, int end) {
    return chunks.append(getRefBases(reference, contig, begin, end));
  }

  private static String getRefBases(FastaReader.Callback.FastaFile reference, String contig,
      int begin, int end) {
    return reference.get(contig, begin - 1, end - 1,
        FastaReader.Callback.FastaFile.Orientation.FORWARD);
  }

  static Optional<String> getSequence(FastaReader.Callback.FastaFile reference,
      String contig, Window window, List<VariantProto> variants) {
    StringBuilder builder = new StringBuilder();
    int homOffset = window.lowerBound();
    for (VariantProto variant : variants) {
      String referenceBases = variant.getReferenceBases();
      int position = variant.getPosition(), next = referenceBases.length() + position;
      if (!Objects.equals(referenceBases, getRefBases(reference, contig, position, next))) {
        return Optional.absent();
      }
      addRefBasesUntil(builder, reference, contig, homOffset, position).append(
          variant.getAlternateBasesList().get(0));
      homOffset = next;
    }
    return Optional.of(
        addRefBasesUntil(builder, reference, contig, homOffset, window.upperBound()).toString());
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
    List<List<VariantProto>>
        falseNegativesQueue = extractVariantQueues(this.falseNegatives, window, location),
        falsePositivesQueue = extractVariantQueues(this.falsePositives, window, location);
    final List<VariantProto> truePositives =
        extractRangeAndFilter(this.truePositives, window, location);
    if (falseNegativesQueue.isEmpty()
        || falsePositivesQueue.isEmpty()
        || WINDOW_MAX_OVERLAPPING < falseNegativesQueue.size() * falsePositivesQueue.size()) {
      return Optional.absent();
    }
    for (final List<VariantProto> newTruePositives : falseNegativesQueue) {
      for (final List<VariantProto> removeFalsePositives : falsePositivesQueue) {
        if (Iterables.any(Iterables.concat(newTruePositives, removeFalsePositives), IS_NON_SNP)) {
          return addTruePosToQueue(newTruePositives, truePositives)
              .transform(
                  new Function<List<VariantProto>, Optional<RescuedVariants>>() {
                    @Override public Optional<RescuedVariants> apply(
                        final List<VariantProto> falseNegatives) {
                      return addTruePosToQueue(removeFalsePositives, truePositives)
                          .transform(
                              new Function<List<VariantProto>, Optional<RescuedVariants>>() {
                                @Override public Optional<RescuedVariants> apply(
                                    List<VariantProto> falsePositives) {
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

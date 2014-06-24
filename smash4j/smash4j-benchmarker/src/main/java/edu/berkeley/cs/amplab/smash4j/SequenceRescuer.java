package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.collect.AbstractSequentialIterator;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.fastaparser.FastaReader.Callback.FastaFile;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.Collection;
import java.util.NavigableMap;
import java.util.Objects;

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
        NavigableMap<Integer, VariantProto> truthLocations,
        NavigableMap<Integer, VariantProto> predictedLocations) {
      return new RescuedVariants(truthLocations, predictedLocations);
    }

    private final NavigableMap<Integer, VariantProto> predictedLocations;
    private final NavigableMap<Integer, VariantProto> truthLocations;

    private RescuedVariants(
        NavigableMap<Integer, VariantProto> truthLocations,
        NavigableMap<Integer, VariantProto> predictedLocations) {
      this.truthLocations = truthLocations;
      this.predictedLocations = predictedLocations;
    }

    public NavigableMap<Integer, VariantProto> newTruePositives() {
      return predictedLocations;
    }

    public NavigableMap<Integer, VariantProto> removeFalsePositives() {
      return truthLocations;
    }
  }

  private static class Window {

    static class Factory {

      static class Builder {

        private static final Ordering<VariantProto>
            LOWEST_START = Ordering.natural().onResultOf(GET_START),
            HIGHEST_END = Ordering.natural().onResultOf(GET_END).reverse();

        private static final int WINDOW_VARIANT_LOOKBACK_SIZE = 50;

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

        private static Function<Window, Window> windowEnlarger(
            final NavigableMap<Integer, VariantProto> variants) {
          return new Function<Window, Window>() {

                @Override public Window apply(Window window) {
                  final int
                      lowerBound = window.lowerBound(),
                      upperBound = window.upperBound();
                  return Window.create(
                      getChoppedVariant(lowerBound, LOWEST_START)
                          .transform(GET_START)
                          .or(lowerBound),
                      getChoppedVariant(upperBound - 1, HIGHEST_END)
                          .transform(
                              new Function<VariantProto, Integer>() {
                                @Override public Integer apply(VariantProto variant) {
                                  int end = GET_END.apply(variant);
                                  return upperBound == end ? end + 1 : end;
                                }
                              })
                          .or(upperBound));
                }

                private Optional<VariantProto> getChoppedVariant(
                    final int location, Ordering<VariantProto> ordering) {
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
              };
        }

        private NavigableMap<Integer, VariantProto> falseNegatives;
        private NavigableMap<Integer, VariantProto> falsePositives;
        private int size;
        private NavigableMap<Integer, VariantProto> truePositives;

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

      private static final int WINDOW_SIZE_LIMIT = 5000;

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
              return variant.getPosition() + variant.getReferenceBases().length();
            }
          };

  public static Builder builder() {
    return new Builder();
  }

  private static Predicate<VariantProto> overlaps(final int location) {
    return new Predicate<VariantProto>() {
      @Override public boolean apply(VariantProto variant) {
        return GET_START.apply(variant) <= location
            && location < GET_END.apply(variant);
      }
    };
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

  public Optional<RescuedVariants> tryRescue(VariantProto variant) {
    throw new UnsupportedOperationException();
  }
}

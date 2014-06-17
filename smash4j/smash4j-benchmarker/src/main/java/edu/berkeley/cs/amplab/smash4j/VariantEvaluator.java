package edu.berkeley.cs.amplab.smash4j;

import static com.google.api.client.repackaged.com.google.common.base.Preconditions.checkNotNull;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Preconditions;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Ordering;
import com.google.common.collect.PeekingIterator;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class VariantEvaluator {

  public static class ContigStats {

    static class Builder {

      private static final
          Function<ImmutableMultiset.Builder<Classification>, Multiset<Classification>>
          build = buildFunction();

      private static <X> Function<ImmutableMultiset.Builder<X>, Multiset<X>> buildFunction() {
        return new Function<ImmutableMultiset.Builder<X>, Multiset<X>>() {
              @Override public Multiset<X> apply(ImmutableMultiset.Builder<X> multiset) {
                return multiset.build();
              }
            };
      }

      private final String contig;
      private final int maxSvBreakpointDistance, rescueWindowSize;
      private final FastaReader.Callback.FastaFile reference;
      private final Map<VariantType, ImmutableMultiset.Builder<Classification>>
          stats = new HashMap<>();
      private final ImmutableList.Builder<VariantProto> variants = ImmutableList.builder();

      private Builder(
          String contig,
          FastaReader.Callback.FastaFile reference,
          int maxSvBreakpointDistance,
          int rescueWindowSize) {
        this.contig = contig;
        this.reference = reference;
        this.maxSvBreakpointDistance = maxSvBreakpointDistance;
        this.rescueWindowSize = rescueWindowSize;
      }

      Builder update(
          long position,
          List<VariantProto> trueVariants,
          List<VariantProto> predictedVariants,
          List<VariantProto> knownFalsePositives) {
        throw new UnsupportedOperationException();
      }

      ContigStats build() {
        return new ContigStats(
            contig,
            Collections.unmodifiableMap(Maps.transformValues(stats, build)),
            variants.build());
      }
    }

    public enum Classification {
      FALSE_NEGATIVE,
      FALSE_POSITIVE,
      KNOWN_FALSE_POSITIVE,
      RESCUED,
      TRUE_POSITIVE
    }

    public enum VariantType {
      INDEL_DELETIONS,
      INDEL_INSERTIONS,
      INDEL_INVERSIONS,
      INDEL_OTHER,
      SNP,
      SV_DELETIONS,
      SV_INSERTIONS,
      SV_OTHER
    }

    static final Function<ContigStats, String> CONTIG =
        new Function<ContigStats, String>() {
          @Override public String apply(ContigStats stats) {
            return stats.contig();
          }
        };

    static Builder builder(
        String contig,
        FastaReader.Callback.FastaFile reference,
        int maxSvBreakpointDistance,
        int rescueWindowSize) {
      return new Builder(contig, reference, maxSvBreakpointDistance, rescueWindowSize);
    }

    private final String contig;
    private final Map<VariantType, Multiset<Classification>> stats;
    private final List<VariantProto> variants;

    private ContigStats(
        String contig,
        Map<VariantType, Multiset<Classification>> stats,
        List<VariantProto> variants) {
      this.contig = contig;
      this.stats = stats;
      this.variants = variants;
    }

    public String contig() {
      return contig;
    }

    public Map<VariantType, Multiset<Classification>> stats() {
      return stats;
    }

    public List<VariantProto> variants() {
      return variants;
    }
  }

  private static class Partitioner<X> {

    static <X> Partitioner<X> fromComparator(Comparator<? super X> comparator) {
      return new Partitioner<>(comparator);
    }

    private final Comparator<? super X> comparator;

    private Partitioner(Comparator<? super X> comparator) {
      this.comparator = comparator;
    }

    FluentIterable<List<X>> partition(final Iterable<? extends X> stream) {
      return new FluentIterable<List<X>>() {
            @Override public Iterator<List<X>> iterator() {
              return new AbstractIterator<List<X>>() {

                    private final PeekingIterator<? extends X>
                        iterator = Iterators.peekingIterator(stream.iterator());

                    @Override protected List<X> computeNext() {
                      if (iterator.hasNext()) {
                        X first = iterator.next();
                        ImmutableList.Builder<X> next = ImmutableList.builder();
                        for (
                            next.add(first);
                            iterator.hasNext() && 0 == comparator.compare(first, iterator.peek());
                            next.add(iterator.next()));
                        return next.build();
                      }
                      return endOfData();
                    }
                  };
            }
          };
    }
  }

  private static class TaggedVariant {

    enum Tag {

      KNOWN_FALSE_POSITIVE,
      PREDICTED,
      TRUE;

      final Function<VariantProto, TaggedVariant> tagger =
          new Function<VariantProto, TaggedVariant>() {
            @Override public TaggedVariant apply(VariantProto variant) {
              return new TaggedVariant(variant, Tag.this);
            }
          };
    }

    static final Function<TaggedVariant, Tag> TAG =
        new Function<TaggedVariant, Tag>() {
          @Override public Tag apply(TaggedVariant taggedVariant) {
            return taggedVariant.tag;
          }
        };

    static final Function<TaggedVariant, VariantProto> VARIANT =
        new Function<TaggedVariant, VariantProto>() {
          @Override public VariantProto apply(TaggedVariant taggedVariant) {
            return taggedVariant.variant;
          }
        };

    private final Tag tag;
    private final VariantProto variant;

    private TaggedVariant(VariantProto variant, Tag tag) {
      this.variant = variant;
      this.tag = tag;
    }
  }

  private static final Function<TaggedVariant, String> GET_CONTIG = Functions.compose(
      new Function<VariantProto, String>() {
        @Override public String apply(VariantProto variant) {
          return variant.getContig();
        }
      },
      TaggedVariant.VARIANT);

  private static final Function<TaggedVariant, Long> GET_POSITION = Functions.compose(
      new Function<VariantProto, Long>() {
        @Override public Long apply(VariantProto variant) {
          return variant.getPosition();
        }
      },
      TaggedVariant.VARIANT);

  private static final Ordering<TaggedVariant>
      CONTIG_ORDERING = Ordering.natural().onResultOf(GET_CONTIG),
      POSITION_ORDERING = Ordering.natural().onResultOf(GET_POSITION),
      CONTIG_POSITION_ORDERING = CONTIG_ORDERING.compound(POSITION_ORDERING);

  private static final Partitioner<TaggedVariant>
      CONTIG_PARTITIONER = Partitioner.fromComparator(CONTIG_ORDERING),
      POSITION_PARTITIONER = Partitioner.fromComparator(POSITION_ORDERING);

  private static <X> X getFirst(
      List<TaggedVariant> taggedVariants, Function<? super TaggedVariant, ? extends X> accessor) {
    return accessor.apply(Preconditions.checkNotNull(Iterables.getFirst(taggedVariants, null)));
  }

  private final Function<List<TaggedVariant>, ContigStats> createStats =
      new Function<List<TaggedVariant>, ContigStats>() {
        @Override public ContigStats apply(List<TaggedVariant> entireContig) {
          ContigStats.Builder stats = ContigStats.builder(
              getFirst(entireContig, GET_CONTIG),
              reference,
              maxSvBreakpointDistance,
              rescueWindowSize);
          for (List<TaggedVariant> variantsAtPos : POSITION_PARTITIONER.partition(entireContig)) {
            ListMultimap<TaggedVariant.Tag, VariantProto> multimap = Multimaps.transformValues(
                Multimaps.index(variantsAtPos, TaggedVariant.TAG),
                TaggedVariant.VARIANT);
            stats.update(
                getFirst(variantsAtPos, GET_POSITION),
                multimap.get(TaggedVariant.Tag.TRUE),
                multimap.get(TaggedVariant.Tag.PREDICTED),
                multimap.get(TaggedVariant.Tag.KNOWN_FALSE_POSITIVE));
          }
          return stats.build();
        }
      };

  private Iterable<VariantProto> knownFalsePositives;
  private int maxSvBreakpointDistance = 100;
  private Iterable<VariantProto> predictedVariants;
  private FastaReader.Callback.FastaFile reference;
  private int rescueWindowSize = 50;
  private Iterable<VariantProto> trueVariants;

  public VariantEvaluator setKnownFalsePositives(Iterable<VariantProto> knownFalsePositives) {
    this.knownFalsePositives = Preconditions.checkNotNull(knownFalsePositives);
    return this;
  }

  public VariantEvaluator setMaxSvBreakpointDistance(int maxSvBreakpointDistance) {
    this.maxSvBreakpointDistance = maxSvBreakpointDistance;
    return this;
  }

  public VariantEvaluator setPredictedVariants(Iterable<VariantProto> predictedVariants) {
    this.predictedVariants = Preconditions.checkNotNull(predictedVariants);
    return this;
  }

  public VariantEvaluator setReference(FastaReader.Callback.FastaFile reference) {
    this.reference = Preconditions.checkNotNull(reference);
    return this;
  }

  public VariantEvaluator setRescueWindowSize(int rescueWindowSize) {
    this.rescueWindowSize = rescueWindowSize;
    return this;
  }

  public VariantEvaluator setTrueVariants(Iterable<VariantProto> trueVariants) {
    this.trueVariants = Preconditions.checkNotNull(trueVariants);
    return this;
  }

  public Map<String, ContigStats> evaluate() {
    checkNotNull(knownFalsePositives, "VariantEvaluator#setKnownFalsePositives was never called");
    checkNotNull(predictedVariants, "VariantEvaluator#setPredictedVariants was never called");
    checkNotNull(reference, "VariantEvaluator#setReference was never called");
    checkNotNull(trueVariants, "VariantEvaluator#setTrueVariants was never called");
    return CONTIG_PARTITIONER
        .partition(Iterables.mergeSorted(
            Arrays.asList(
                FluentIterable.from(trueVariants)
                    .transform(TaggedVariant.Tag.TRUE.tagger),
                FluentIterable.from(predictedVariants)
                    .transform(TaggedVariant.Tag.PREDICTED.tagger),
                FluentIterable.from(knownFalsePositives)
                    .transform(TaggedVariant.Tag.KNOWN_FALSE_POSITIVE.tagger)),
            CONTIG_POSITION_ORDERING))
        .transform(createStats)
        .uniqueIndex(ContigStats.CONTIG);
  }
}

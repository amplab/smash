package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Normalizer {

  private static final Logger LOGGER = Logger.getLogger(Normalizer.class.getName());

  private class PositionCollisionHandler implements Iterable<PositionCollisionHandler.Contig> {

    class Contig implements Iterable<Contig.Position> {

      class Position implements Iterable<Variant> {

        private final Map.Entry<Integer, Collection<Variant>> entry;

        private Position(Map.Entry<Integer, Collection<Variant>> entry) {
          this.entry = entry;
        }

        @Override public Iterator<Variant> iterator() {
          Collection<Variant> value = entry.getValue();
          if (1 == value.size()) {
            return value.iterator();
          }
          List<Variant>
              normed = new ArrayList<>(),
              notnormed = new ArrayList<>();
          OUTER: for (Variant variant : value) {
            for (String key : variant.info().keySet()) {
              if (NORM_INFO_TAG.equals(key)) {
                normed.add(variant);
                continue OUTER;
              }
            }
            notnormed.add(variant);
          }
          Variant basevar;
          if (!notnormed.isEmpty()) {
            basevar = Preconditions.checkNotNull(Iterables.getFirst(notnormed, null));
            for (Variant variant : Iterables.skip(notnormed, 1)) {
              String chrom = variant.contig();
              long pos = variant.position();
              LOGGER.warning(String.format(
                  "Variant already exists on %s at %d; discarding variant %s %d %s/%s",
                  chrom,
                  pos,
                  chrom,
                  pos,
                  variant.referenceBases(),
                  Joiner.on(",").join(variant.alternateBases())));
            }
          } else {
            basevar = Preconditions.checkNotNull(Iterables.getFirst(normed, null));
            normed = ImmutableList.copyOf(Iterables.skip(normed, 1));
          }
          List<Variant> finalVars = new ArrayList<>();
          for (Variant variant : normed) {
            try {
              finalVars.add(basevar);
              basevar = shiftUntilNotOverlapping(basevar, variant);
            } catch (AssertionError e) {
              Variant first = Preconditions.checkNotNull(Iterables.getFirst(normed, null));
              LOGGER.warning(
                  String.format("failed denorm at %s %d", first.contig(), first.position()));
            }
          }
          finalVars.add(basevar);
          return finalVars.iterator();
        }

        private Variant shiftUntilNotOverlapping(
            Variant varOne, final Variant varTwo) {
          final String varTwoContig = varTwo.contig();
          int onePos = varOne.position() - 1,
              twoPos = varTwo.position() - 1;
          String twoRefAllele = varTwo.referenceBases();
          List<String> twoAltAlleles = varTwo.alternateBases();
          for (int varOneLastPos = onePos + varOne.referenceBases().length();
              twoPos < varOneLastPos;) {
            final int twoPosCopy = ++twoPos;
            final String twoRefAlleleCopy = twoRefAllele;
            twoRefAllele = slide(varTwoContig, twoPos, twoRefAllele, twoRefAllele);
            twoAltAlleles = FluentIterable.from(twoAltAlleles)
                .transform(
                    new Function<String, String>() {
                      @Override public String apply(String slidingAllele) {
                        return slide(
                            varTwoContig, twoPosCopy, twoRefAlleleCopy, slidingAllele);
                      }
                    })
                .toList();
            assert sameFirstAndLastBase(twoRefAllele, twoAltAlleles);
          }
          return varTwo.toBuilder()
              .setPosition(twoPos + 1)
              .setReferenceBases(twoRefAllele)
              .setAlternateBases(twoAltAlleles)
              .build();
        }

        private String slide(String contig, int pos, String refAllele, String slidingAllele) {
          int rightEndpoint = pos + refAllele.length() - 1;
          return slidingAllele.substring(1)
              + fastaFile.get(
                  contig,
                  rightEndpoint,
                  rightEndpoint + 1,
                  FastaReader.Callback.FastaFile.Orientation.FORWARD);
        }

        private boolean sameFirstAndLastBase(String refAllele, List<String> altAlleles) {
          Collection<String> allAlleles =
              ImmutableList.copyOf(Iterables.concat(Collections.singleton(refAllele), altAlleles));
          boolean b = true;
          for (SameBaseTester tester : SameBaseTester.values()) {
            b &= tester.sameBase(allAlleles);
          }
          return b;
        }
      }

      private final Map.Entry<String, Map<Integer, Collection<Variant>>> entry;

      final Function<Map.Entry<Integer, Collection<Variant>>, Position> newPosition =
          new Function<Map.Entry<Integer, Collection<Variant>>, Position>() {
            @Override public Position apply(Map.Entry<Integer, Collection<Variant>> entry) {
              return new Position(entry);
            }
          };

      private Contig(Map.Entry<String, Map<Integer, Collection<Variant>>> entry) {
        this.entry = entry;
      }

      @Override public Iterator<Position> iterator() {
        return Iterators.transform(entry.getValue().entrySet().iterator(), newPosition);
      }
    }

    private final Map<String, Map<Integer, Collection<Variant>>> cache;

    final Function<Map.Entry<String, Map<Integer, Collection<Variant>>>, Contig> newContig =
        new Function<Map.Entry<String, Map<Integer, Collection<Variant>>>, Contig>() {
          @Override public Contig apply(
              Map.Entry<String, Map<Integer, Collection<Variant>>> entry) {
            return new Contig(entry);
          }
        };

    private PositionCollisionHandler(Map<String, Map<Integer, Collection<Variant>>> cache) {
      this.cache = cache;
    }

    @Override public Iterator<Contig> iterator() {
      return Iterators.transform(cache.entrySet().iterator(), newContig);
    }
  }

  private enum SameBaseTester {

    FIRST_BASE {
      @Override int index(String string) {
        return 0;
      }
    },

    LAST_BASE {
      @Override int index(String string) {
        return string.length() - 1;
      }
    };

    abstract int index(String string);

    final boolean sameBase(Iterable<String> strings) {
      Character base = null;
      for (String string : strings) {
        if (string.isEmpty()) {
          return false;
        }
        char c = string.charAt(index(string));
        if (null == base) {
          base = c;
        } else if (!base.equals(c)) {
          return false;
        }
      }
      return null != base;
    }
  }

  private static class VariantFilter implements Predicate<Variant> {

    private static final Pattern SYMBOLIC_ALLELE = Pattern.compile("<.+?>");

    static VariantFilter create(int maxIndelSize) {
      return new VariantFilter(maxIndelSize);
    }

    private final int maxIndelSize;

    private VariantFilter(int maxIndelSize) {
      this.maxIndelSize = maxIndelSize;
    }

    @Override public boolean apply(Variant variant) {
      if (maxIndelSize < variant.referenceBases().length()) {
        return true;
      }
      for (String alternateBases : variant.alternateBases()) {
        if (maxIndelSize < alternateBases.length()
            || SYMBOLIC_ALLELE.matcher(alternateBases).matches()) {
          return true;
        }
      }
      return Optional.fromNullable(Iterables.getFirst(Genotype.getGenotypes(variant), null))
          .transform(fromPredicate(Predicates.not(Predicates.or(
              Predicates.equalTo(Genotype.NO_CALL),
              Predicates.equalTo(Genotype.HOM_REF)))))
          .or(Boolean.FALSE);
    }

    private static <X> Function<X, Boolean> fromPredicate(final Predicate<? super X> predicate) {
      return new Function<X, Boolean>() {
            @Override public Boolean apply(X object) {
              return predicate.apply(object) ? Boolean.TRUE : Boolean.FALSE;
            }
          };
    }
  }

  static final String NORM_INFO_TAG = "OP";

  private static final Function<String, String> TO_UPPER_CASE =
      new Function<String, String>() {
        @Override public String apply(String string) {
          return string.toUpperCase();
        }
      };

  public static Normalizer cleanOnly(int maxIndelSize, FastaReader.Callback.FastaFile fastaFile) {
    return new Normalizer(maxIndelSize, fastaFile, true);
  }

  public static Normalizer create(int maxIndelSize, FastaReader.Callback.FastaFile fastaFile) {
    return new Normalizer(maxIndelSize, fastaFile, false);
  }

  private static Function<String, String> redundancyChopper(
      final String ref, final List<String> alts) {
    return new Function<String, String>() {

          private final Pattern pattern = createPattern(
              ImmutableList.copyOf(Iterables.concat(Collections.singletonList(ref), alts)));

          @Override public String apply(String string) {
            Matcher matcher = pattern.matcher(string);
            return matcher.matches() ? matcher.group(1) : string;
          }

          private Pattern createPattern(List<String> strings) {
            StringBuilder builder = new StringBuilder();
            OUTER: for (int i = 1; true; ++i) {
              Character character = null;
              for (String string : strings) {
                int j = string.length() - i;
                if (j < 1) {
                  break OUTER;
                }
                char c = string.charAt(j);
                if (null == character) {
                  character = c;
                } else if (!character.equals(c)) {
                  break OUTER;
                }
              }
              if (null == character) {
                break;
              }
              builder.append(character);
            }
            return Pattern.compile(
                String.format("(.*?)%s", Pattern.quote(builder.reverse().toString())));
          }
        };
  }

  private final boolean cleanOnly;
  private final FastaReader.Callback.FastaFile fastaFile;
  private final int maxIndelSize;

  private final Function<Variant, Variant> normalize =
      new Function<Variant, Variant>() {

        @Override public Variant apply(Variant variant) {
          String ref = TO_UPPER_CASE.apply(variant.referenceBases());
          List<String> alts = FluentIterable.from(variant.alternateBases())
              .transform(TO_UPPER_CASE)
              .toList();
          if (cleanOnly) {
            return variant.toBuilder()
                .setReferenceBases(ref)
                .setAlternateBases(alts)
                .build();
          }
          Function<String, String> chopper = redundancyChopper(ref, alts);
          ref = chopper.apply(ref);
          alts = FluentIterable.from(alts)
              .transform(chopper)
              .toList();
          int originalPosition = variant.position(), pos = originalPosition - 1;
          for (String contig = variant.contig(); SameBaseTester.LAST_BASE.sameBase(
              Iterables.concat(Collections.singletonList(ref), alts));) {
            final String prevBase = TO_UPPER_CASE.apply(fastaFile.get(
                contig,
                --pos,
                pos + 1,
                FastaReader.Callback.FastaFile.Orientation.FORWARD));
            Function<String, String> slider =
                new Function<String, String>() {
                  @Override public String apply(String string) {
                    return new StringBuilder(prevBase)
                        .append(string.substring(0, string.length() - 1))
                        .toString();
                  }
                };
            ref = slider.apply(ref);
            alts = FluentIterable.from(alts).transform(slider).toList();
          }
          int newPosition = pos + 1;
          Variant.Builder builder = variant.toBuilder()
              .setPosition(newPosition)
              .setReferenceBases(ref)
              .setAlternateBases(alts);
          if (newPosition < originalPosition) {
            builder.setInfo(
                ImmutableListMultimap.<String, String>builder()
                    .putAll(variant.info())
                    .put(NORM_INFO_TAG, String.valueOf(originalPosition))
                    .build());
          }
          return builder.build();
        }
      };

  private Normalizer(
      int maxIndelSize, FastaReader.Callback.FastaFile fastaFile, boolean cleanOnly) {
    this.maxIndelSize = maxIndelSize;
    this.fastaFile = fastaFile;
    this.cleanOnly = cleanOnly;
  }

  private PositionCollisionHandler createPositionCollisionHandler(Iterable<Variant> variants) {
    Map<String, Map<Integer, Collection<Variant>>> cache = new TreeMap<>();
    for (Variant variant : variants) {
      String contig = variant.contig();
      int position = variant.position();
      Map<Integer, Collection<Variant>> map = cache.get(contig);
      if (null == map) {
        cache.put(contig, map = new TreeMap<>());
      }
      Collection<Variant> list = map.get(position);
      if (null == list) {
        map.put(position, list = new ArrayList<>());
      }
      list.add(variant);
    }
    return new PositionCollisionHandler(cache);
  }

  public FluentIterable<Variant> normalize(Iterable<Variant> variants) {
    FluentIterable<Variant> intermediate = FluentIterable.from(variants)
        .filter(VariantFilter.create(maxIndelSize))
        .transform(normalize);
    return cleanOnly
        ? intermediate
        : FluentIterable.from(
            Iterables.concat(Iterables.concat(createPositionCollisionHandler(intermediate))));
  }
}

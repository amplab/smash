package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

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

  private class PositionCollisionFilter implements Iterable<PositionCollisionFilter.Contig> {

    class Contig implements Iterable<Contig.Position> {

      class Position implements Iterable<VariantProto> {

        private final Map.Entry<Long, Collection<VariantProto>> entry;

        private Position(Map.Entry<Long, Collection<VariantProto>> entry) {
          this.entry = entry;
        }

        @Override public Iterator<VariantProto> iterator() {
          Collection<VariantProto> value = entry.getValue();
          if (1 == value.size()) {
            return value.iterator();
          }
          List<VariantProto>
              normed = new ArrayList<>(),
              notnormed = new ArrayList<>();
          OUTER: for (VariantProto variant : value) {
            for (VariantProto.Multimap.Entry entry : variant.getInfo().getEntryList()) {
              if (NORM_INFO_TAG.equals(entry.getKey())) {
                normed.add(variant);
                continue OUTER;
              }
            }
            notnormed.add(variant);
          }
          VariantProto basevar;
          if (!notnormed.isEmpty()) {
            basevar = Preconditions.checkNotNull(Iterables.getFirst(notnormed, null));
            for (VariantProto variant : Iterables.skip(notnormed, 1)) {
              String chrom = variant.getContig();
              long pos = variant.getPosition();
              LOGGER.warning(String.format(
                  "Variant already exists on %s at %d; discarding variant %s %d %s/%s",
                  chrom,
                  pos,
                  chrom,
                  pos,
                  variant.getReferenceBases(),
                  Joiner.on(",").join(variant.getAlternateBasesList())));
            }
          } else {
            basevar = Preconditions.checkNotNull(Iterables.getFirst(normed, null));
            normed = ImmutableList.copyOf(Iterables.skip(normed, 1));
          }
          List<VariantProto> finalVars = new ArrayList<>();
          for (VariantProto variant : normed) {
            try {
              finalVars.add(basevar);
              basevar = shiftUntilNotOverlapping(basevar, variant);
            } catch (AssertionError e) {
              VariantProto first = Preconditions.checkNotNull(Iterables.getFirst(normed, null));
              LOGGER.warning(
                  String.format("failed denorm at %s %d", first.getContig(), first.getPosition()));
            }
          }
          finalVars.add(basevar);
          return finalVars.iterator();
        }

        private VariantProto shiftUntilNotOverlapping(
            VariantProto varOne, final VariantProto varTwo) {
          final String varTwoContig = varTwo.getContig();
          int
              onePos = (int) varOne.getPosition() - 1,
              twoPos = (int) varTwo.getPosition() - 1;
          String twoRefAllele = varTwo.getReferenceBases();
          List<String> twoAltAlleles = varTwo.getAlternateBasesList();
          for (int varOneLastPos = onePos + varOne.getReferenceBases().length();
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
              .clearAlternateBases()
              .addAllAlternateBases(twoAltAlleles)
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

      private final Map.Entry<String, Map<Long, Collection<VariantProto>>> entry;

      final Function<Map.Entry<Long, Collection<VariantProto>>, Position> newPosition =
          new Function<Map.Entry<Long, Collection<VariantProto>>, Position>() {
            @Override public Position apply(Map.Entry<Long, Collection<VariantProto>> entry) {
              return new Position(entry);
            }
          };

      private Contig(Map.Entry<String, Map<Long, Collection<VariantProto>>> entry) {
        this.entry = entry;
      }

      @Override public Iterator<Position> iterator() {
        return Iterators.transform(entry.getValue().entrySet().iterator(), newPosition);
      }
    }

    private final Map<String, Map<Long, Collection<VariantProto>>> cache;

    final Function<Map.Entry<String, Map<Long, Collection<VariantProto>>>, Contig> newContig =
        new Function<Map.Entry<String, Map<Long, Collection<VariantProto>>>, Contig>() {
          @Override public Contig apply(
              Map.Entry<String, Map<Long, Collection<VariantProto>>> entry) {
            return new Contig(entry);
          }
        };

    private PositionCollisionFilter(Map<String, Map<Long, Collection<VariantProto>>> cache) {
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

  private static class VariantFilter implements Predicate<VariantProto> {

    private static final Function<VariantProto.Multimap, Boolean> NOT_HOMREF =
        new Function<VariantProto.Multimap, Boolean>() {
          @Override public Boolean apply(VariantProto.Multimap call) {
            for (VariantProto.Multimap.Entry entry : call.getEntryList()) {
              if ("GT".equals(entry.getKey())) {
                switch (Iterables.getOnlyElement(entry.getValueList())) {
                  case ".":
                  case "0/0":
                  case "0|0":
                    return false;
                  default:
                    return true;
                }
              }
            }
            throw new IllegalStateException("Missing \"GT\" INFO value");
          }
        };

    private static final Pattern SYMBOLIC_ALLELE = Pattern.compile("<.+?>");

    static VariantFilter create(int maxIndelSize) {
      return new VariantFilter(maxIndelSize);
    }

    private final int maxIndelSize;

    private VariantFilter(int maxIndelSize) {
      this.maxIndelSize = maxIndelSize;
    }

    @Override public boolean apply(VariantProto variant) {
      if (maxIndelSize < variant.getReferenceBases().length()) {
        return true;
      }
      for (String alternateBases : variant.getAlternateBasesList()) {
        if (maxIndelSize < alternateBases.length()
            || SYMBOLIC_ALLELE.matcher(alternateBases).matches()) {
          return true;
        }
      }
      return Optional.fromNullable(Iterables.getFirst(variant.getCallList(), null))
          .transform(NOT_HOMREF)
          .or(Boolean.FALSE);
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

  private final Function<VariantProto, VariantProto> normalize =
      new Function<VariantProto, VariantProto>() {

        @Override public VariantProto apply(VariantProto variant) {
          String ref = TO_UPPER_CASE.apply(variant.getReferenceBases());
          List<String> alts = FluentIterable.from(variant.getAlternateBasesList())
              .transform(TO_UPPER_CASE)
              .toList();
          if (cleanOnly) {
            return variant.toBuilder()
                .setReferenceBases(ref)
                .clearAlternateBases()
                .addAllAlternateBases(alts)
                .build();
          }
          Function<String, String> chopper = redundancyChopper(ref, alts);
          ref = chopper.apply(ref);
          alts = FluentIterable.from(alts)
              .transform(chopper)
              .toList();
          int originalPosition = (int) variant.getPosition(), pos = originalPosition - 1;
          for (String contig = variant.getContig(); SameBaseTester.LAST_BASE.sameBase(
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
          VariantProto.Builder builder = variant.toBuilder()
              .setPosition(newPosition)
              .setReferenceBases(ref)
              .clearAlternateBases()
              .addAllAlternateBases(alts);
          if (newPosition < originalPosition) {
            builder.getInfoBuilder().addEntry(VariantProto.Multimap.Entry.newBuilder()
                .setKey(NORM_INFO_TAG)
                .addValue(String.valueOf(originalPosition)));
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

  private PositionCollisionFilter createPositionCollisionFilter(Iterable<VariantProto> variants) {
    Map<String, Map<Long, Collection<VariantProto>>> cache = new TreeMap<>();
    for (VariantProto variant : variants) {
      String contig = variant.getContig();
      long position = variant.getPosition();
      Map<Long, Collection<VariantProto>> map = cache.get(contig);
      if (null == map) {
        cache.put(contig, map = new TreeMap<>());
      }
      Collection<VariantProto> list = map.get(position);
      if (null == list) {
        map.put(position, list = new ArrayList<>());
      }
      list.add(variant);
    }
    return new PositionCollisionFilter(cache);
  }

  public FluentIterable<VariantProto> normalize(Iterable<VariantProto> variants) {
    return FluentIterable.from(Iterables.concat(Iterables.concat(createPositionCollisionFilter(
        FluentIterable.from(variants)
            .filter(VariantFilter.create(maxIndelSize))
            .transform(normalize)))));
  }
}

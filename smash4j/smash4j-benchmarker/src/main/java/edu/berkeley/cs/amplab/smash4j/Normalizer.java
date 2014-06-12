package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Normalizer {

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

  private static final Function<String, String> TO_UPPER_CASE =
      new Function<String, String>() {
        @Override public String apply(String string) {
          return string.toUpperCase();
        }
      };

  static final String NORM_INFO_TAG = "OP";

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
          for (
              String contig = variant.getContig();
              sameLastBase(Iterables.concat(Collections.singletonList(ref), alts));) {
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

        private boolean sameLastBase(Iterable<String> strings) {
          Character lastBase = null;
          for (String string : strings) {
            if (string.isEmpty()) {
              return false;
            }
            char c = string.charAt(string.length() - 1);
            if (null == lastBase) {
              lastBase = c;
            } else if (!lastBase.equals(c)) {
              return false;
            }
          }
          return null != lastBase;
        }
      };

  private Normalizer(
      int maxIndelSize, FastaReader.Callback.FastaFile fastaFile, boolean cleanOnly) {
    this.maxIndelSize = maxIndelSize;
    this.fastaFile = fastaFile;
    this.cleanOnly = cleanOnly;
  }

  public FluentIterable<VariantProto> normalize(Iterable<VariantProto> variants) {
    return FluentIterable.from(variants)
        .filter(VariantFilter.create(maxIndelSize))
        .transform(normalize);
  }
}

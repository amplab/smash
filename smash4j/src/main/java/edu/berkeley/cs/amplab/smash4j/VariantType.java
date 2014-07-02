package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;

public enum VariantType {

  INDEL_DELETION(false, false),
  INDEL_INSERTION(false, false),
  INDEL_INVERSION(false, false),
  INDEL_OTHER(false, false),
  SNP(true, false),
  SV_DELETION(false, true),
  SV_INSERTION(false, true),
  SV_OTHER(false, true);

  private static String getFirstAlt(Variant variant) {
    return Iterables.getOnlyElement(variant.alternateBases());
  }

  public static Function<Variant, VariantType> getType(final int maxIndelSize) {
    return new Function<Variant, VariantType>() {
          @Override public VariantType apply(Variant variant) {
            return isSnp(variant)
                ? SNP
                : isStructuralVariant(variant, maxIndelSize)
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
        };
  }

  private static boolean hasSingleAlt(Variant variant) {
    return 1 == variant.alternateBases().size();
  }

  private static boolean isDeletion(Variant variant) {
    return getFirstAlt(variant).length() < variant.referenceBases().length();
  }

  private static boolean isInsertion(Variant variant) {
    return variant.referenceBases().length() < getFirstAlt(variant).length();
  }

  private static boolean isInversion(Variant variant) {
    String ref = variant.referenceBases(), alt = getFirstAlt(variant);
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

  private static boolean isSnp(Variant variant) {
    if (1 == variant.referenceBases().length()) {
      for (String alternateBases : variant.alternateBases()) {
        if (1 != alternateBases.length()) {
          return false;
        }
      }
      return true;
    }
    return false;
  }

  private static boolean isStructuralVariant(Variant variant, int maxIndelSize) {
    if (variant.referenceBases().length() <= maxIndelSize) {
      for (String alternateBases : variant.alternateBases()) {
        if (maxIndelSize < alternateBases.length()) {
          return true;
        }
      }
      for (String key : variant.info().keySet()) {
        if ("SVTYPE".equals(key)) {
          return true;
        }
      }
    }
    return false;
  }

  private final boolean isSnp;
  private final boolean isStructuralVariant;

  private VariantType(boolean isSnp, boolean isStructuralVariant) {
    this.isSnp = isSnp;
    this.isStructuralVariant = isStructuralVariant;
  }

  public boolean isSnp() {
    return isSnp;
  }

  public boolean isStructuralVariant() {
    return isStructuralVariant;
  }
}

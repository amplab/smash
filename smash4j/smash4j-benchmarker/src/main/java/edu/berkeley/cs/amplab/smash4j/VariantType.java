package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

public enum VariantType {

  INDEL_DELETION(false, false),
  INDEL_INSERTION(false, false),
  INDEL_INVERSION(false, false),
  INDEL_OTHER(false, false),
  SNP(true, false),
  SV_DELETION(false, true),
  SV_INSERTION(false, true),
  SV_OTHER(false, true);

  private static String getFirstAlt(VariantProto variant) {
    return Iterables.getOnlyElement(variant.getAlternateBasesList());
  }

  public static Function<VariantProto, VariantType> getType(final int maxIndelSize) {
    return new Function<VariantProto, VariantType>() {
          @Override public VariantType apply(VariantProto variant) {
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

  private static boolean isStructuralVariant(VariantProto variant, int maxIndelSize) {
    if (variant.getReferenceBases().length() <= maxIndelSize) {
      for (String alternateBases : variant.getAlternateBasesList()) {
        if (maxIndelSize < alternateBases.length()) {
          return true;
        }
      }
      for (VariantProto.Multimap.Entry entry : variant.getInfo().getEntryList()) {
        if ("SVTYPE".equals(entry.getKey())) {
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

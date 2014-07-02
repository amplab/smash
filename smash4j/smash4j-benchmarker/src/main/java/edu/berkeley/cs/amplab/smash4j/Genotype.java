package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.collect.FluentIterable;

import java.util.Iterator;
import java.util.List;

public enum Genotype {

  HET,
  HOM_REF,
  HOM_VAR,
  NO_CALL;

  public static final Function<Variant, List<Genotype>> GET_GENOTYPES =
      new Function<Variant, List<Genotype>>() {
        @Override public List<Genotype> apply(Variant variant) {
          return getGenotypes(variant);
        }
      };

  private static final Function<Variant.Call, Genotype> GET_GENOTYPE =
      new Function<Variant.Call, Genotype>() {
        @Override public Genotype apply(Variant.Call call) {
          Integer i = null;
          for (Iterator<Integer> iterator = call.genotype().iterator(); iterator.hasNext();) {
            Integer j = iterator.next();
            if (null == i) {
              i = j;
            } else if (!i.equals(j)) {
              return HET;
            }
          }
          return null != i ? i.equals(0) ? HOM_REF : HOM_VAR : NO_CALL;
        }
      };

  public static List<Genotype> getGenotypes(Variant variant) {
    return FluentIterable.from(variant.calls()).transform(GET_GENOTYPE).toList();
  }
}

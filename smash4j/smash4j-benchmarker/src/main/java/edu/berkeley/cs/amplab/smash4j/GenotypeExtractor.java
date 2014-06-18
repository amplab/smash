package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Iterables;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.List;

public enum GenotypeExtractor {

  INSTANCE;

  private static final Function<VariantProto.Multimap, String> GET_GENOTYPE =
      new Function<VariantProto.Multimap, String>() {

        private final Predicate<VariantProto.Multimap.Entry> isGt = Predicates.compose(
            Predicates.equalTo("GT"),
            new Function<VariantProto.Multimap.Entry, String>() {
              @Override public String apply(VariantProto.Multimap.Entry entry) {
                return entry.getKey();
              }
            });

        private final Function<VariantProto.Multimap.Entry, List<String>> getValueList =
            new Function<VariantProto.Multimap.Entry, List<String>>() {
              @Override public List<String> apply(VariantProto.Multimap.Entry entry) {
                return entry.getValueList();
              }
            };

        @Override public String apply(VariantProto.Multimap call) {
          return Iterables.getOnlyElement(FluentIterable.from(call.getEntryList()).filter(isGt)
              .transformAndConcat(getValueList));
        }
      };

  public Optional<String> getGenotype(VariantProto variant) {
    return Optional.fromNullable(Iterables.getFirst(variant.getCallList(), null))
        .transform(GET_GENOTYPE);
  }
}

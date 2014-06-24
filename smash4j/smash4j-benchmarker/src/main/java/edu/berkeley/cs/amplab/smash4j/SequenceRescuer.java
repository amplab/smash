package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.fastaparser.FastaReader.Callback.FastaFile;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.util.NavigableMap;

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
          Preconditions.checkNotNull(rescueWindowSize));
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

  public static Builder builder() {
    return new Builder();
  }

  private final String contig;
  private final NavigableMap<Integer, VariantProto> falseNegatives;
  private final NavigableMap<Integer, VariantProto> falsePositives;
  private final FastaReader.Callback.FastaFile reference;
  private final int rescueWindowSize;
  private final NavigableMap<Integer, VariantProto> truePositives;

  private SequenceRescuer(
      String contig,
      NavigableMap<Integer, VariantProto> truePositives,
      NavigableMap<Integer, VariantProto> falsePositives,
      NavigableMap<Integer, VariantProto> falseNegatives,
      FastaFile reference,
      int rescueWindowSize) {
    this.contig = contig;
    this.truePositives = truePositives;
    this.falsePositives = falsePositives;
    this.falseNegatives = falseNegatives;
    this.reference = reference;
    this.rescueWindowSize = rescueWindowSize;
  }

  public Optional<RescuedVariants> rescue(VariantProto variant) {
    throw new UnsupportedOperationException();
  }
}

package edu.berkeley.cs.amplab.vcfparser;

import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public final class VcfRecord {

  public static final class Builder {

    private Optional<List<String>> alt = Optional.absent();
    private Optional<String> chrom = Optional.absent();
    private Optional<List<String>> filters = Optional.absent();
    private Optional<List<String>> format = Optional.absent();
    private Optional<List<String>> ids = Optional.absent();
    private Optional<Map<String, List<String>>> info = Optional.absent();
    private Optional<Integer> pos = Optional.absent();
    private Optional<Double> qual = Optional.absent();
    private Optional<String> ref = Optional.absent();
    private final ImmutableList.Builder<List<String>> samples = ImmutableList.builder();

    private Builder() {}

    private Builder(
        String chrom,
        Integer pos,
        List<String> ids,
        String ref,
        List<String> alt,
        Double qual,
        List<String> filters,
        Map<String, List<String>> info,
        List<String> format,
        List<List<String>> samples) {
      setChrom(chrom);
      setPos(pos);
      setIds(ids);
      setRef(ref);
      setAlt(alt);
      setQual(qual);
      setFilters(filters);
      setInfo(info);
      setFormat(format);
      for (List<String> sample : samples) {
        addSample(sample);
      }
    }

    public Builder setChrom(String chrom) {
      this.chrom = Optional.fromNullable(chrom);
      return this;
    }

    public Builder setPos(int pos) {
      this.pos = Optional.fromNullable(pos);
      return this;
    }

    public Builder setIds(List<String> ids) {
      this.ids = Optional.fromNullable(ids);
      return this;
    }

    public Builder setRef(String ref) {
      this.ref = Optional.fromNullable(ref);
      return this;
    }

    public Builder setAlt(List<String> alt) {
      this.alt = Optional.fromNullable(alt);
      return this;
    }

    public Builder setQual(double qual) {
      this.qual = Optional.fromNullable(qual);
      return this;
    }

    public Builder setFilters(List<String> filters) {
      this.filters = Optional.fromNullable(filters);
      return this;
    }

    public Builder setInfo(Map<String, List<String>> info) {
      this.info = Optional.fromNullable(info);
      return this;
    }

    public Builder setFormat(List<String> format) {
      this.format = Optional.fromNullable(format);
      return this;
    }

    public Builder addSample(List<String> sample) {
      this.samples.add(sample);
      return this;
    }

    public VcfRecord build() {
      List<String> format = this.format.orNull();
      List<List<String>> samples = this.samples.build();
      Preconditions.checkState(
          null == format == samples.isEmpty(),
          "If there are no samples, format can't be set. If there are any samples, format must be set");
      return new VcfRecord(
          chrom.orNull(),
          pos.orNull(),
          ids.orNull(),
          ref.orNull(),
          alt.orNull(),
          qual.orNull(),
          filters.orNull(),
          info.orNull(),
          format,
          samples);
    }
  }

  public static Builder builder() {
    return new Builder();
  }

  private final List<String> alt;
  private final String chrom;
  private final List<String> filters;
  private final List<String> format;
  private final List<String> ids;
  private final Map<String, List<String>> info;
  private final Integer pos;
  private final Double qual;
  private final String ref;
  private final List<List<String>> samples;

  private VcfRecord(
      String chrom,
      Integer pos,
      List<String> ids,
      String ref,
      List<String> alt,
      Double qual,
      List<String> filters,
      Map<String, List<String>> info,
      List<String> format,
      List<List<String>> samples) {
    this.chrom = chrom;
    this.pos = pos;
    this.ids = ids;
    this.ref = ref;
    this.alt = alt;
    this.qual = qual;
    this.filters = filters;
    this.info = info;
    this.format = format;
    this.samples = samples;
  }

  @Override
  public boolean equals(Object obj) {
    boolean same = this == obj;
    if (!same && null != obj && VcfRecord.class == obj.getClass()) {
      VcfRecord rhs = (VcfRecord) obj;
      return Objects.equals(chrom(), rhs.chrom())
          && Objects.equals(pos(), rhs.pos())
          && Objects.equals(ids(), rhs.ids())
          && Objects.equals(ref(), rhs.ref())
          && Objects.equals(alt(), rhs.alt())
          && Objects.equals(qual(), rhs.qual())
          && Objects.equals(filters(), rhs.filters())
          && Objects.equals(info(), rhs.info())
          && Objects.equals(format(), rhs.format())
          && Objects.equals(samples(), rhs.samples());
    }
    return same;
  }

  public String chrom() {
    return chrom;
  }

  public int pos() {
    return pos;
  }

  public List<String> ids() {
    return ids;
  }

  public String ref() {
    return ref;
  }

  public List<String> alt() {
    return alt;
  }

  public double qual() {
    return qual;
  }

  public List<String> filters() {
    return filters;
  }

  public Map<String, List<String>> info() {
    return info;
  }

  public List<String> format() {
    return format;
  }

  public List<List<String>> samples() {
    return samples;
  }

  @Override
  public int hashCode() {
    return Objects.hash(
        chrom(),
        pos(),
        ids(),
        ref(),
        alt(),
        qual(),
        filters(),
        info(),
        format(),
        samples());
  }

  public Builder toBuilder() {
    return new Builder(
        chrom(),
        pos(),
        ids(),
        ref(),
        alt(),
        qual(),
        filters(),
        info(),
        format(),
        samples());
  }

  @Override
  public String toString() {
    Optional<List<String>> format = Optional.fromNullable(format());
    return String.format(
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s",
        toString(chrom()),
        toString(pos()),
        toString(ids(), ';'),
        toString(ref()),
        toString(alt(), ','),
        toString(qual()),
        toString(filters(), ';'),
        toString(info()),
        format.isPresent() ? String.format("\t%s", Joiner.on(':').join(format.get())) : "",
        toString(samples()));
  }

  private static String toString(Object object) {
    return null == object ? "." : Objects.toString(object);
  }

  private static String toString(List<String> list, char separator) {
    return null == list ? "." : Joiner.on(separator).join(list);
  }

  private static String toString(Map<String, List<String>> map) {
    if (null != map) {
      StringBuilder builder = new StringBuilder();
      for (
          Iterator<Map.Entry<String, List<String>>> iterator = map.entrySet().iterator();
          iterator.hasNext();
          builder.append(iterator.hasNext() ? ";" : "")) {
        Map.Entry<String, List<String>> entry = iterator.next();
        List<String> value = entry.getValue();
        builder.append(String.format(
            "%s%s",
            entry.getKey(),
            value.isEmpty() ? "" : String.format("=%s", toString(value, ','))));
      }
      return builder.toString();
    }
    return ".";
  }

  private static String toString(List<List<String>> samples) {
    if (!samples.isEmpty()) {
      StringBuilder builder = new StringBuilder("\t");
      for (
          Iterator<List<String>> iterator = samples.iterator();
          iterator.hasNext();
          builder.append(Joiner.on(':').join(iterator.next()))
              .append(iterator.hasNext() ? "\t" : ""));
      return builder.toString();
    }
    return "";
  }
}

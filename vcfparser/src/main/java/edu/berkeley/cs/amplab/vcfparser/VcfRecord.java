package edu.berkeley.cs.amplab.vcfparser;

import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;

import java.util.List;
import java.util.Objects;

public final class VcfRecord {

  private final String chrom;
  private final Integer pos;
  private final String id;
  private final String ref;
  private final String alt;
  private final Integer qual;
  private final String filter;
  private final String info;
  private final String format;
  private final List<String> samples;

  private VcfRecord(
      String chrom,
      Integer pos,
      String id,
      String ref,
      String alt,
      Integer qual,
      String filter,
      String info,
      String format,
      List<String> samples) {
    this.chrom = chrom;
    this.pos = pos;
    this.id = id;
    this.ref = ref;
    this.alt = alt;
    this.qual = qual;
    this.filter = filter;
    this.info = info;
    this.format = format;
    this.samples = samples;
  }

  public static Builder builder() {
    return new Builder();
  }

  public Builder toBuilder() {
    return new Builder(
        chrom(),
        pos(),
        id(),
        ref(),
        alt(),
        qual(),
        filter(),
        info(),
        format(),
        samples());
  }

  public String chrom() {
    return chrom;
  }

  public int pos() {
    return pos;
  }

  public String id() {
    return id;
  }

  public String ref() {
    return ref;
  }

  public String alt() {
    return alt;
  }

  public int qual() {
    return qual;
  }

  public String filter() {
    return filter;
  }

  public String info() {
    return info;
  }

  public String format() {
    return format;
  }

  public List<String> samples() {
    return samples;
  }

  @Override
  public int hashCode() {
    return Objects.hash(
        chrom(),
        pos(),
        id(),
        ref(),
        alt(),
        qual(),
        filter(),
        info(),
        format(),
        samples());
  }

  @Override
  public boolean equals(Object obj) {
    boolean same = this == obj;
    if (!same && null != obj && VcfRecord.class == obj.getClass()) {
      VcfRecord rhs = (VcfRecord) obj;
      return Objects.equals(chrom(), rhs.chrom())
          && Objects.equals(pos(), rhs.pos())
          && Objects.equals(id(), rhs.id())
          && Objects.equals(ref(), rhs.ref())
          && Objects.equals(alt(), rhs.alt())
          && Objects.equals(qual(), rhs.qual())
          && Objects.equals(filter(), rhs.filter())
          && Objects.equals(info(), rhs.info())
          && Objects.equals(format(), rhs.format())
          && Objects.equals(samples(), rhs.samples());
    }
    return same;
  }

  @Override
  public String toString() {
    Optional<String> format = Optional.fromNullable(format());
    List<String> samples = samples();
    return String.format(
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s",
        toString(chrom()),
        toString(pos()),
        toString(id()),
        toString(ref()),
        toString(alt()),
        toString(qual()),
        toString(filter()),
        toString(info()),
        format.isPresent() ? String.format("\t%s", format.get()) : "",
        samples.isEmpty() ? "" : String.format("\t%s", Joiner.on('\t').join(samples)));
  }

  private static String toString(Object object) {
    return null == object ? "." : Objects.toString(object);
  }

  public static final class Builder {

    private final ImmutableList.Builder<String> samples = ImmutableList.builder();
    private Optional<String> chrom = Optional.absent();
    private Optional<Integer> pos = Optional.absent();
    private Optional<String> id = Optional.absent();
    private Optional<String> ref = Optional.absent();
    private Optional<String> alt = Optional.absent();
    private Optional<Integer> qual = Optional.absent();
    private Optional<String> filter = Optional.absent();
    private Optional<String> info = Optional.absent();
    private Optional<String> format = Optional.absent();

    private Builder() {}

    private Builder(
        String chrom,
        Integer pos,
        String id,
        String ref,
        String alt,
        Integer qual,
        String filter,
        String info,
        String format,
        List<String> samples) {
      setChrom(chrom);
      setPos(pos);
      setId(id);
      setRef(ref);
      setAlt(alt);
      setQual(qual);
      setFilter(filter);
      setInfo(info);
      setFormat(format);
      for (String sample : samples) {
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

    public Builder setId(String id) {
      this.id = Optional.fromNullable(id);
      return this;
    }

    public Builder setRef(String ref) {
      this.ref = Optional.fromNullable(ref);
      return this;
    }

    public Builder setAlt(String alt) {
      this.alt = Optional.fromNullable(alt);
      return this;
    }

    public Builder setQual(int qual) {
      this.qual = Optional.fromNullable(qual);
      return this;
    }

    public Builder setFilter(String filter) {
      this.filter = Optional.fromNullable(filter);
      return this;
    }

    public Builder setInfo(String info) {
      this.info = Optional.fromNullable(info);
      return this;
    }

    public Builder setFormat(String format) {
      this.format = Optional.fromNullable(format);
      return this;
    }

    public Builder addSample(String sample) {
      this.samples.add(sample);
      return this;
    }

    public VcfRecord build() {
      return new VcfRecord(
          chrom.orNull(),
          pos.orNull(),
          id.orNull(),
          ref.orNull(),
          alt.orNull(),
          qual.orNull(),
          filter.orNull(),
          info.orNull(),
          format.orNull(),
          samples.build());
    }
  }
}

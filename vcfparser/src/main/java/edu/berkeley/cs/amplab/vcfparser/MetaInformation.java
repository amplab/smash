package edu.berkeley.cs.amplab.vcfparser;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;

import java.net.URL;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * The meta-information records in a VCF file.
 */
public final class MetaInformation {

  private final Set<MetaInfoLine<?>> lines;

  private MetaInformation(Set<MetaInfoLine<?>> lines) {
    this.lines = lines;
  }

  public static Builder builder(FileFormat.Format format) {
    return new Builder(Collections.singleton(FileFormat.create(format)));
  }

  public String toString() {
    StringBuilder builder = new StringBuilder();
    for (MetaInfoLine<?> line : lines) {
      builder.append(String.format("%s%n", line));
    }
    return builder.toString();
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(lines);
  }

  @Override
  public boolean equals(Object obj) {
    return this == obj
        || null != obj
        && MetaInformation.class == obj.getClass()
        && Objects.equals(lines, ((MetaInformation) obj).lines);
  }

  public Builder toBuilder() {
    return new Builder(lines);
  }

  public FileFormat.Format format() {
    for (MetaInfoLine<?> line : lines) {
      if (line instanceof FileFormat) {
        return ((FileFormat) line).format();
      }
    }
    throw new IllegalStateException("Format line not found");
  }

  /**
   * Return the {@code MetaInformation.Info} records contained in the meta-information.
   */
  public Set<Info> infos() {
    return get(Info.class);
  }

  /**
   * Return the {@code MetaInformation.Format} records contained in the meta-information.
   */
  public Set<Format> formats() {
    return get(Format.class);
  }

  /**
   * Return the {@code MetaInformation.Filter} records contained in the meta-information.
   */
  public Set<Filter> filters() {
    return get(Filter.class);
  }

  /**
   * Return the {@code MetaInformation.Alt} records contained in the meta-information.
   */
  public Set<Alt> alts() {
    return get(Alt.class);
  }

  /**
   * Return the {@code MetaInformation.Assembly} records contained in the meta-information.
   */
  public Set<Assembly> assemblies() {
    return get(Assembly.class);
  }

  /**
   * Return the {@code MetaInformation.Contig} records contained in the meta-information.
   */
  public Set<Contig> contigs() {
    return get(Contig.class);
  }

  /**
   * Return the {@code MetaInformation.Sample} records contained in the meta-information.
   */
  public Set<Sample> samples() {
    return get(Sample.class);
  }

  /**
   * Return the {@code MetaInformation.Pedigree} records contained in the meta-information.
   */
  public Set<Pedigree> pedigrees() {
    return get(Pedigree.class);
  }

  /**
   * Return the {@code MetaInformation.PedigreeDB} records contained in the meta-information.
   */
  public Set<PedigreeDB> pedigreeDBs() {
    return get(PedigreeDB.class);
  }

  /**
   * Return the {@code MetaInformation.UnparsedMetaInfoLine} records contained in the meta-information.
   */
  public Set<UnparsedMetaInfoLine> unparsedLines() {
    return get(UnparsedMetaInfoLine.class);
  }

  private <X> Set<X> get(Class<X> type) {
    ImmutableSet.Builder<X> infos = ImmutableSet.builder();
    for (MetaInfoLine<?> line : lines) {
      if (type.isInstance(line)) {
        infos.add(type.cast(line));
      }
    }
    return infos.build();
  }

  private abstract static class EscapingString {

    private static final Pattern ESCAPE_PATTERN = Pattern.compile("\\\\|\"");
    private final String value;

    private EscapingString(String value) {
      this.value = value;
    }

    public static EscapingString escaping(Object value) {
      return
          new EscapingString(value.toString()) {
            @Override
            public String toString() {
              return String.format("\"%s\"", ESCAPE_PATTERN.matcher(value()).replaceAll("\\\\$0"));
            }
          };
    }

    public static EscapingString nonEscaping(Object value) {
      return
          new EscapingString(value.toString()) {
            @Override
            public String toString() {
              return value();
            }
          };
    }

    String value() {
      return value;
    }

    @Override
    public final int hashCode() {
      return Objects.hashCode(value());
    }

    @Override
    public final boolean equals(Object obj) {
      return this == obj
          || null != obj
          && getClass() == obj.getClass()
          && Objects.equals(value(), ((EscapingString) obj).value());
    }

    @Override
    public abstract String toString();
  }

  private abstract static class MetaInfoLine<V> {

    private final String type;
    private final V value;

    private MetaInfoLine(String type, V value) {
      this.type = type;
      this.value = value;
    }

    String lineType() {
      return type;
    }

    V value() {
      return value;
    }

    @Override
    public final int hashCode() {
      return Objects.hash(type, value());
    }

    @Override
    public final boolean equals(Object obj) {
      boolean same = this == obj;
      if (!same && null != obj && getClass() == obj.getClass()) {
        MetaInfoLine<?> rhs = (MetaInfoLine<?>) obj;
        return Objects.equals(type, rhs.type)
            && Objects.equals(value(), rhs.value());
      }
      return same;
    }

    @Override
    public final String toString() {
      return String.format("##%s=%s", type, toStringImpl(value()));
    }

    abstract String toStringImpl(V value);
  }

  private abstract static class SimpleMetaInfoLine<V> extends MetaInfoLine<V> {

    private SimpleMetaInfoLine(String type, V value) {
      super(type, value);
    }

    @Override
    public final String toStringImpl(V value) {
      return Objects.toString(value);
    }
  }

  private abstract static class CompoundMetaInfoLine extends MetaInfoLine<Map<String, EscapingString>> {

    private CompoundMetaInfoLine(String type, Map<String, EscapingString> value) {
      super(type, value);
    }

    @Override
    public final String toStringImpl(Map<String, EscapingString> value) {
      StringBuilder builder = new StringBuilder("<");
      for (
          Iterator<Map.Entry<String, EscapingString>> iterator = value.entrySet().iterator();
          iterator.hasNext();) {
        Map.Entry<String, EscapingString> next = iterator.next();
        builder.append(next.getKey())
            .append('=')
            .append(next.getValue())
            .append(iterator.hasNext() ? "," : "");
      }
      return builder.append('>').toString();
    }

    String get(String key) {
      return value().get(key).value();
    }

    Map<String, String> extraFields(Set<String> requiredFields) {
      ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      for (Map.Entry<String, EscapingString> entry : value().entrySet()) {
        String key = entry.getKey();
        if (!requiredFields.contains(key)) {
          extraFields.put(key, entry.getValue().value());
        }
      }
      return extraFields.build();
    }
  }

  /**
   * The "fileformat" metainfo record.
   */
  public static final class FileFormat extends SimpleMetaInfoLine<FileFormat.Format> {

    private FileFormat(Format format) {
      super("fileformat", format);
    }

    public static FileFormat create(Format format) {
      return new FileFormat(format);
    }

    public Format format() {
      return value();
    }

    public enum Format {

      V4_0,
      V4_1,
      V4_2;

      public static Format parse(String value) {
        switch (value) {
          case "VCFv4.0":
            return V4_0;
          case "VCFv4.1":
            return V4_1;
          case "VCFv4.2":
            return V4_2;
          default:
            throw new IllegalArgumentException(String.format("Failed to parse file format: \"%s\"", value));
        }
      }

      @Override
      public String toString() {
        return String.format("VCFv%s", name().substring(1).replace('_', '.'));
      }
    }
  }

  /**
   * The "assembly" metainfo record.
   */
  public static final class Assembly extends SimpleMetaInfoLine<URL> {

    private Assembly(URL url) {
      super("assembly", url);
    }

    public static Assembly create(URL url) {
      return new Assembly(url);
    }

    public URL url() {
      return value();
    }
  }

  /**
   * The metainfo record for record types that the parser doesn't know how to parse into rich objects.
   */
  public static final class UnparsedMetaInfoLine extends SimpleMetaInfoLine<String> {

    private UnparsedMetaInfoLine(String type, String value) {
      super(type, value);
    }

    public static UnparsedMetaInfoLine create(String type, String value) {
      return new UnparsedMetaInfoLine(type, value);
    }

    public String type() {
      return lineType();
    }

    public String value() {
      return super.value();
    }
  }

  /**
   * The type of the @{code Number} field in {@code MetaInformation.Info} and {@code MetaInformation.Format} records.
   */
  public static final class Number {

    public static final Number
        A = create('A'),
        R = create('R'),
        G = create('G'),
        UNKNOWN = create('.');

    private final String value;

    private Number(String value) {
      this.value = value;
    }

    public static Number create(int i) {
      Preconditions.checkArgument(0 <= i, "Number must be nonnegative");
      return new Number(String.valueOf(i));
    }

    public static Number create(String value) {
      switch (value) {
        case "A":
          return A;
        case "R":
          return R;
        case "G":
          return G;
        case ".":
          return UNKNOWN;
        default:
          return create(Integer.parseInt(value));
      }
    }

    private static Number create(char c) {
      return new Number(String.valueOf(c));
    }

    @Override
    public int hashCode() {
      return Objects.hashCode(toString());
    }

    @Override
    public boolean equals(Object obj) {
      return this == obj
          || null != obj
          && Number.class == obj.getClass()
          && Objects.equals(toString(), obj.toString());
    }

    @Override
    public String toString() {
      return value;
    }
  }

  /**
   * The "INFO" metainfo record.
   */
  public static final class Info extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID", "Number", "Type", "Description");

    private Info(Map<String, EscapingString> value) {
      super("INFO", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public Number number() {
      return Number.create(get("Number"));
    }

    public Type type() {
      return Type.parse(get("Type"));
    }

    public String description() {
      return get("Description");
    }

    public Map<String, String> extraFields() {
      return extraFields(REQUIRED_FIELDS);
    }

    public Builder toBuilder() {
      return new Builder(id(), number(), type(), description(), extraFields());
    }

    public enum Type {

      INTEGER,
      FLOAT,
      FLAG,
      CHARACTER,
      STRING;

      public static Type parse(String value) {
        return valueOf(value.toUpperCase());
      }

      @Override
      public String toString() {
        String name = name();
        return name.substring(0, 1) + name.substring(1).toLowerCase();
      }
    }

    /**
     * The builder class for {@code MetaInformation.Info} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      private Optional<String> id = Optional.absent();
      private Optional<Number> number = Optional.absent();
      private Optional<Type> type = Optional.absent();
      private Optional<String> description = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          Number number,
          Type type,
          String description,
          Map<String, String> extraFields) {
        this.id = Optional.of(id);
        this.number = Optional.of(number);
        this.type = Optional.of(type);
        this.description = Optional.of(description);
        this.extraFields.putAll(extraFields);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder setNumber(Number number) {
        this.number = Optional.of(number);
        return this;
      }

      public Builder setType(Type type) {
        this.type = Optional.of(type);
        return this;
      }

      public Builder setDescription(String description) {
        this.description = Optional.of(description);
        return this;
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Info build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()))
            .put("Number", EscapingString.nonEscaping(number.get()))
            .put("Type", EscapingString.nonEscaping(type.get()))
            .put("Description", EscapingString.escaping(description.get()));
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.escaping(entry.getValue()));
        }
        return new Info(value.build());
      }
    }
  }

  /**
   * The "FILTER" metainfo record.
   */
  public static final class Filter extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID", "Description");

    private Filter(Map<String, EscapingString> value) {
      super("FILTER", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public String description() {
      return get("Description");
    }

    public Map<String, String> extraFields() {
      return extraFields(REQUIRED_FIELDS);
    }

    public Builder toBuilder() {
      return new Builder(id(), description(), extraFields());
    }

    /**
     * The builder class for {@code MetaInformation.Filter} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      private Optional<String> id = Optional.absent();
      private Optional<String> description = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          String description,
          Map<String, String> extraFields) {
        this.id = Optional.of(id);
        this.description = Optional.of(description);
        this.extraFields.putAll(extraFields);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder setDescription(String description) {
        this.description = Optional.of(description);
        return this;
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Filter build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()))
            .put("Description", EscapingString.escaping(description.get()));
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.escaping(entry.getValue()));
        }
        return new Filter(value.build());
      }
    }
  }

  /**
   * The "FORMAT" metainfo record.
   */
  public static final class Format extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID", "Number", "Type", "Description");

    private Format(Map<String, EscapingString> value) {
      super("FORMAT", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public Number number() {
      return Number.create(get("Number"));
    }

    public Type type() {
      return Type.parse(get("Type"));
    }

    public String description() {
      return get("Description");
    }

    public Map<String, String> extraFields() {
      return extraFields(REQUIRED_FIELDS);
    }

    public Builder toBuilder() {
      return new Builder(id(), number(), type(), description(), extraFields());
    }

    public enum Type {

      INTEGER,
      FLOAT,
      CHARACTER,
      STRING;

      public static Type parse(String value) {
        return valueOf(value.toUpperCase());
      }

      @Override
      public String toString() {
        String name = name();
        return name.substring(0, 1) + name.substring(1).toLowerCase();
      }
    }

    /**
     * The builder class for {@code MetaInformation.Format} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      private Optional<String> id = Optional.absent();
      private Optional<Number> number = Optional.absent();
      private Optional<Type> type = Optional.absent();
      private Optional<String> description = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          Number number,
          Type type,
          String description,
          Map<String, String> extraFields) {
        this.id = Optional.of(id);
        this.number = Optional.of(number);
        this.type = Optional.of(type);
        this.description = Optional.of(description);
        this.extraFields.putAll(extraFields);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder setNumber(Number number) {
        this.number = Optional.of(number);
        return this;
      }

      public Builder setType(Type type) {
        this.type = Optional.of(type);
        return this;
      }

      public Builder setDescription(String description) {
        this.description = Optional.of(description);
        return this;
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Format build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()))
            .put("Number", EscapingString.nonEscaping(number.get()))
            .put("Type", EscapingString.nonEscaping(type.get()))
            .put("Description", EscapingString.escaping(description.get()));
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.escaping(entry.getValue()));
        }
        return new Format(value.build());
      }
    }
  }

  /**
   * The "ALT" metainfo record.
   */
  public static final class Alt extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID", "Description");

    private Alt(Map<String, EscapingString> value) {
      super("ALT", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public String description() {
      return get("Description");
    }

    public Map<String, String> extraFields() {
      return extraFields(REQUIRED_FIELDS);
    }

    public Builder toBuilder() {
      return new Builder(id(), description(), extraFields());
    }

    /**
     * The builder class for {@code MetaInformation.Alt} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      private Optional<String> id = Optional.absent();
      private Optional<String> description = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          String description,
          Map<String, String> extraFields) {
        this.id = Optional.of(id);
        this.description = Optional.of(description);
        this.extraFields.putAll(extraFields);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder setDescription(String description) {
        this.description = Optional.of(description);
        return this;
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Alt build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()))
            .put("Description", EscapingString.escaping(description.get()));
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.escaping(entry.getValue()));
        }
        return new Alt(value.build());
      }
    }
  }

  /**
   * The "contig" metainfo record.
   */
  public static final class Contig extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID");

    private Contig(Map<String, EscapingString> value) {
      super("contig", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public Map<String, String> extraFields() {
      return extraFields(REQUIRED_FIELDS);
    }

    public Builder toBuilder() {
      return new Builder(id(), extraFields());
    }

    /**
     * The builder class for {@code MetaInformation.Contig} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();
      private Optional<String> id = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          Map<String, String> extraFields) {
        this.id = Optional.of(id);
        this.extraFields.putAll(extraFields);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Contig build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()));
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.nonEscaping(entry.getValue()));
        }
        return new Contig(value.build());
      }
    }
  }

  /**
   * The "SAMPLE" metainfo record.
   */
  public static final class Sample extends CompoundMetaInfoLine {

    static final Set<String> REQUIRED_FIELDS = ImmutableSet.of("ID", "Genome", "Mixture", "Description");

    private Sample(Map<String, EscapingString> value) {
      super("SAMPLE", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public String id() {
      return get("ID");
    }

    public String genome() {
      return get("Genome");
    }

    public String mixture() {
      return get("Mixture");
    }

    public String description() {
      return get("Description");
    }

    public Builder toBuilder() {
      return new Builder(id(), genome(), mixture(), description());
    }

    /**
     * The builder class for {@code MetaInformation.Sample} objects.
     */
    public static final class Builder {

      private Optional<String> id = Optional.absent();
      private Optional<String> genome = Optional.absent();
      private Optional<String> mixture = Optional.absent();
      private Optional<String> description = Optional.absent();

      private Builder() {}

      private Builder(
          String id,
          String genome,
          String mixture,
          String description) {
        this.id = Optional.of(id);
        this.genome = Optional.of(genome);
        this.mixture = Optional.of(mixture);
        this.description = Optional.of(description);
      }

      public Builder setId(String id) {
        this.id = Optional.of(id);
        return this;
      }

      public Builder setGenome(String genome) {
        this.genome = Optional.of(genome);
        return this;
      }

      public Builder setMixture(String mixture) {
        this.mixture = Optional.of(mixture);
        return this;
      }

      public Builder setDescription(String description) {
        this.description = Optional.of(description);
        return this;
      }

      public Sample build() {
        return new Sample(ImmutableMap.<String, EscapingString>builder()
            .put("ID", EscapingString.nonEscaping(id.get()))
            .put("Genome", EscapingString.nonEscaping(genome.get()))
            .put("Mixture", EscapingString.nonEscaping(mixture.get()))
            .put("Description", EscapingString.nonEscaping(description.get())).build());
      }
    }
  }

  /**
   * The "PEDIGREE" metainfo record.
   */
  public static final class Pedigree extends CompoundMetaInfoLine {

    private Pedigree(Map<String, EscapingString> value) {
      super("PEDIGREE", value);
    }

    public static Builder builder() {
      return new Builder();
    }

    public Map<String, String> extraFields() {
      return extraFields(Collections.<String>emptySet());
    }

    public Builder toBuilder() {
      return new Builder(extraFields());
    }

    /**
     * The builder class for {@code MetaInformation.Pedigree} objects.
     */
    public static final class Builder {

      private final ImmutableMap.Builder<String, String> extraFields = ImmutableMap.builder();

      private Builder() {}

      private Builder(Map<String, String> extraFields) {
        this.extraFields.putAll(extraFields);
      }

      public Builder addExtraField(String key, String value) {
        this.extraFields.put(key, value);
        return this;
      }

      public Pedigree build() {
        ImmutableMap.Builder<String, EscapingString> value = ImmutableMap.<String, EscapingString>builder();
        for (Map.Entry<String, String> entry : extraFields.build().entrySet()) {
          value.put(entry.getKey(), EscapingString.nonEscaping(entry.getValue()));
        }
        return new Pedigree(value.build());
      }
    }
  }

  /**
   * The "pedigreeDB" metainfo record.
   */
  public static final class PedigreeDB extends SimpleMetaInfoLine<URL> {

    private PedigreeDB(URL url) {
      super("pedigreeDB", url);
    }

    public static PedigreeDB create(URL url) {
      return new PedigreeDB(url);
    }

    public URL url() {
      return value();
    }
  }

  /**
   * The builder class for {@code MetaInformation} objects.
   */
  public static final class Builder {

    private final ImmutableSet.Builder<MetaInfoLine<?>> lines = ImmutableSet.builder();

    private Builder(Iterable<? extends MetaInfoLine<?>> lines) {
      this.lines.addAll(lines);
    }

    public Builder addInfo(Info info) {
      return add(info);
    }

    public Builder addInfo(Info.Builder info) {
      return addInfo(info.build());
    }

    public Builder addFilter(Filter filter) {
      return add(filter);
    }

    public Builder addFilter(Filter.Builder filter) {
      return addFilter(filter.build());
    }

    public Builder addFormat(Format format) {
      return add(format);
    }

    public Builder addFormat(Format.Builder format) {
      return addFormat(format.build());
    }

    public Builder addAlt(Alt alt) {
      return add(alt);
    }

    public Builder addAlt(Alt.Builder alt) {
      return addAlt(alt.build());
    }

    public Builder addAssembly(Assembly assembly) {
      return add(assembly);
    }

    public Builder addContig(Contig contig) {
      return add(contig);
    }

    public Builder addContig(Contig.Builder contig) {
      return addContig(contig.build());
    }

    public Builder addSample(Sample sample) {
      return add(sample);
    }

    public Builder addSample(Sample.Builder sample) {
      return addSample(sample.build());
    }

    public Builder addPedigree(Pedigree pedigree) {
      return add(pedigree);
    }

    public Builder addPedigree(Pedigree.Builder pedigree) {
      return addPedigree(pedigree.build());
    }

    public Builder addPedigreeDB(PedigreeDB pedigreeDB) {
      return add(pedigreeDB);
    }

    public Builder addUnparsedLine(UnparsedMetaInfoLine line) {
      return add(line);
    }

    private Builder add(MetaInfoLine<?> line) {
      lines.add(line);
      return this;
    }

    public MetaInformation build() {
      return new MetaInformation(lines.build());
    }
  }
}

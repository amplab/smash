package edu.berkeley.cs.amplab.vcfparser;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Throwables;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Iterators;
import com.google.common.collect.Maps;
import com.google.common.collect.PeekingIterator;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A reader for VCF data.
 */
public class VcfReader implements Closeable {

  /**
   * A callback type used to receive VCF data from the source it is read from.
   */
  public interface Callback<X> {

    /**
     * Compute some value from the given {@code MetaInformation}, {@code Header}, and {@code FluentIterable<VcfRecord>}
     * while the data is being consumed from the source.
     */
    X readVcf(MetaInformation metaInformation, Header header, FluentIterable<VcfRecord> records);
  }

  /**
   * Static factory method that creates a new {@code VcfReader} that reads from the given {@code InputStream} source.
   */
  public static VcfReader from(InputStream in) {
    return from(new InputStreamReader(in));
  }

  /**
   * Static factory method that creates a new {@code VcfReader} that reads from the given {@code Writer} source.
   */
  public static VcfReader from(Reader in) {
    return new VcfReader(in instanceof BufferedReader ? (BufferedReader) in : new BufferedReader(in));
  }

  private final BufferedReader in;

  private VcfReader(BufferedReader in) {
    this.in = in;
  }

  private static final Pattern METAINFO_LINE_PATTERN = Pattern.compile("##(\\p{Alpha}+?)=(.+)");

  /**
   * Read the VCF data and pass the data to the provided {@code Callback}.
   */
  public <X> X read(Callback<X> callback) {
    final PeekingIterator<String> iterator = Iterators.peekingIterator(
        new AbstractIterator<String>() {
          @Override
          protected String computeNext() {
            try {
              String line = in.readLine();
              return null == line ? endOfData() : line;
            } catch (IOException e) {
              throw Throwables.propagate(e);
            }
          }
        }
    );
    return callback.readVcf(
        parseMetaInfo(iterator),
        parseHeader(iterator),
        new FluentIterable<VcfRecord>() {
          @Override
          public Iterator<VcfRecord> iterator() {
            return Iterators.transform(iterator, PARSE_RECORD);
          }
        });
  }

  private static final Pattern HEADER_LINE_PATTERN = Pattern.compile("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(\tFORMAT(?:\t\\p{Alnum}+?)+)?");

  private static final Function<String, VcfRecord> PARSE_RECORD =
      new Function<String, VcfRecord>() {
        @Override
        public VcfRecord apply(String line) {
          return parseRecord(line);
        }
      };

  private static VcfRecord parseRecord(String line) {
    String[] parts = line.split("\t");
    VcfRecord.Builder record = VcfRecord.builder();
    if (isDefined(parts[0])) {
      record.setChrom(parts[0]);
    }
    if (isDefined(parts[1])) {
      record.setPos(Integer.parseInt(parts[1]));
    }
    if (isDefined(parts[2])) {
      record.setId(parts[2]);
    }
    if (isDefined(parts[3])) {
      record.setRef(parts[3]);
    }
    if (isDefined(parts[4])) {
      record.setAlt(parts[4]);
    }
    if (isDefined(parts[5])) {
      record.setQual(Integer.parseInt(parts[5]));
    }
    if (isDefined(parts[6])) {
      record.setFilter(parts[6]);
    }
    if (isDefined(parts[7])) {
      record.setInfo(parts[7]);
    }
    if (isDefined(parts[8])) {
      record.setFormat(parts[8]);
    }
    for (int i = 9; i < parts.length; ++i) {
      record.addSample(parts[i]);
    }
    return record.build();
  }

  private static boolean isDefined(String value) {
    return !(null == value || ".".equals(value));
  }

  private static Header parseHeader(Iterator<String> iterator) {
    if (iterator.hasNext()) {
      String next = iterator.next();
      Matcher matcher = HEADER_LINE_PATTERN.matcher(next);
      if (matcher.matches()) {
        String[] parts = Optional.fromNullable(matcher.group(1)).or("").split("\t");
        Header.Builder header = Header.builder();
        for (int i = 2; i < parts.length; ++i) {
          header.addSampleId(parts[i]);
        }
        return header.build();
      }
      throw new IllegalStateException(String.format("Failed to parse header line: \"%s\"", next));
    }
    throw new IllegalStateException("Empty file");
  }

  private static MetaInformation parseMetaInfo(PeekingIterator<String> iterator) {
    MetaInformation.Builder metainfo = MetaInformation.builder(getFormat(iterator));
    while (iterator.hasNext() && iterator.peek().startsWith("##")) {
      String next = iterator.next();
      Matcher matcher = METAINFO_LINE_PATTERN.matcher(next);
      if (!matcher.matches()) {
        throw new IllegalStateException(String.format("Failed to parse metainfo line: \"%s\"", next));
      }
      String type = matcher.group(1);
      String value = matcher.group(2);
      switch (type) {
        case "INFO":
          metainfo.addInfo(info(value));
          break;
        case "FILTER":
          metainfo.addFilter(filter(value));
          break;
        case "FORMAT":
          metainfo.addFormat(format(value));
          break;
        case "ALT":
          metainfo.addAlt(alt(value));
          break;
        case "assembly":
          metainfo.addAssembly(assembly(value));
          break;
        case "contig":
          metainfo.addContig(contig(value));
          break;
        case "SAMPLE":
          metainfo.addSample(sample(value));
          break;
        case "PEDIGREE":
          metainfo.addPedigree(pedigree(value));
          break;
        case "pedigreeDB":
          metainfo.addPedigreeDB(pedigreeDB(value));
          break;
        default:
          metainfo.addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create(type, value));
      }
    }
    return metainfo.build();
  }

  private static MetaInformation.FileFormat.Format getFormat(Iterator<String> iterator) {
    if (iterator.hasNext()) {
      String next = iterator.next();
      Matcher matcher = METAINFO_LINE_PATTERN.matcher(next);
      if (matcher.matches()) {
        if ("fileformat".equals(matcher.group(1))) {
          return MetaInformation.FileFormat.Format.parse(matcher.group(2));
        }
        throw new IllegalStateException("fileformat must be the first record in the file");
      }
      throw new IllegalStateException(String.format("Failed to parse metainfo line: \"%s\"", next));
    }
    throw new IllegalStateException("Empty file");
  }

  private static final Pattern KEY_VALUE_PAIR_PATTERN = Pattern.compile("(\\p{Alnum}+?)=(\"(?:\\\\\"|[^\"])*?\"|[^\"][^,]*?)(,|$)");

  private static MetaInformation.Info info(String value) {
    Map<String, String> map = parseCompoundValue(value);
    MetaInformation.Info.Builder builder = MetaInformation.Info.builder()
        .setId(map.get("ID"))
        .setNumber(MetaInformation.Number.create(map.get("Number")))
        .setType(MetaInformation.Info.Type.parse(map.get("Type")))
        .setDescription(map.get("Description"));
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Info.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    for (Map.Entry<String, String> entry : extraFields.entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.Filter filter(String value) {
    Map<String, String> map = parseCompoundValue(value);
    MetaInformation.Filter.Builder builder = MetaInformation.Filter.builder()
        .setId(map.get("ID"))
        .setDescription(map.get("Description"));
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Filter.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    for (Map.Entry<String, String> entry : extraFields.entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.Format format(String value) {
    Map<String, String> map = parseCompoundValue(value);
    MetaInformation.Format.Builder builder = MetaInformation.Format.builder()
        .setId(map.get("ID"))
        .setNumber(MetaInformation.Number.create(map.get("Number")))
        .setType(MetaInformation.Format.Type.parse(map.get("Type")))
        .setDescription(map.get("Description"));
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Format.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    for (Map.Entry<String, String> entry : extraFields.entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.Alt alt(String value) {
    Map<String, String> map = parseCompoundValue(value);
    MetaInformation.Alt.Builder builder = MetaInformation.Alt.builder()
        .setId(map.get("ID"))
        .setDescription(map.get("Description"));
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Alt.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    for (Map.Entry<String, String> entry : extraFields.entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.Assembly assembly(String value) {
    try {
      return MetaInformation.Assembly.create(new URL(value));
    } catch (MalformedURLException e) {
      throw Throwables.propagate(e);
    }
  }

  private static MetaInformation.Contig contig(String value) {
    Map<String, String> map = parseCompoundValue(value, false);
    MetaInformation.Contig.Builder builder = MetaInformation.Contig.builder()
        .setId(map.get("ID"));
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Contig.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    for (Map.Entry<String, String> entry : extraFields.entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.Sample sample(String value) {
    Map<String, String> map = parseCompoundValue(value);
    return MetaInformation.Sample.builder()
        .setId(map.get("ID"))
        .setGenome(map.get("Genome"))
        .setMixture(map.get("Mixture"))
        .setDescription(map.get("Description"))
        .build();
  }

  private static MetaInformation.Pedigree pedigree(String value) {
    MetaInformation.Pedigree.Builder builder = MetaInformation.Pedigree.builder();
    for (Map.Entry<String, String> entry : parseCompoundValue(value).entrySet()) {
      builder.addExtraField(entry.getKey(), entry.getValue());
    }
    return builder.build();
  }

  private static MetaInformation.PedigreeDB pedigreeDB(String value) {
    try {
      return MetaInformation.PedigreeDB.create(new URL(value));
    } catch (MalformedURLException e) {
      throw Throwables.propagate(e);
    }
  }

  private static Map<String, String> parseCompoundValue(String value) {
    return parseCompoundValue(value, true);
  }

  private static Map<String, String> parseCompoundValue(String value, boolean unescape) {
    if (value.startsWith("<") && value.endsWith(">")) {
      ImmutableMap.Builder<String, String> map = ImmutableMap.builder();
      for (Matcher matcher = KEY_VALUE_PAIR_PATTERN.matcher(value.substring(1, value.length() - 1)); matcher.find();) {
        String group2 = matcher.group(2);
        map.put(matcher.group(1), unescape ? unescape(group2) : group2);
        if (matcher.group(3).isEmpty()) {
          break;
        }
      }
      return map.build();
    }
    throw new IllegalStateException(String.format("Malformed value: \"%s\"", value));
  }

  private static final Pattern UNESCAPE_PATTERN = Pattern.compile("\\\\(\\\\|\")");

  private static String unescape(String value) {
    return value.startsWith("\"") && value.endsWith("\"")
        ? UNESCAPE_PATTERN.matcher(value.substring(1, value.length() - 1)).replaceAll("$1")
        : value;
  }

  @Override
  public void close() throws IOException {
    in.close();
  }
}

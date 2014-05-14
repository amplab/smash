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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A reader for VCF data.
 */
public class VcfReader {

  /**
   * A callback used to receive VCF data from the source it is read from.
   */
  public interface Callback<X> {

    /**
     * Compute a value of type {@code X>} from the given metainfo, header, and records.
     * The value computed will be returned from {@link #read}.
     */
    X readVcf(MetaInformation metaInformation, Header header, FluentIterable<VcfRecord> records);
  }

  private static final Pattern
      HEADER_LINE_PATTERN = Pattern.compile("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(\tFORMAT(?:\t\\p{Alnum}+?)+)?"),
      KEY_VALUE_PAIR_PATTERN = Pattern.compile("(\\p{Alnum}+?)=(\"(?:\\\\\"|[^\"])*?\"|[^\"][^,]*?)(,|$)"),
      METAINFO_LINE_PATTERN = Pattern.compile("##(\\p{Alpha}+?)=(.+)"),
      UNESCAPE_PATTERN = Pattern.compile("\\\\(\\\\|\")");

  private static final Function<String, VcfRecord>
      PARSE_RECORD =
          new Function<String, VcfRecord>() {
            @Override public VcfRecord apply(String line) {
              return parseRecord(line);
            }
          };

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

  private static VcfRecord parseRecord(String line) {
    String
        fields[] = line.split("\t"),
        chrom = fields[0],
        pos = fields[1],
        id = fields[2],
        ref = fields[3],
        alt = fields[4],
        qual = fields[5],
        filter = fields[6],
        info = fields[7],
        format = fields[8];
    VcfRecord.Builder record = VcfRecord.builder();
    if (isDefined(chrom)) {
      record.setChrom(chrom);
    }
    if (isDefined(pos)) {
      record.setPos(Integer.parseInt(pos));
    }
    if (isDefined(id)) {
      record.setId(id);
    }
    if (isDefined(ref)) {
      record.setRef(ref);
    }
    if (isDefined(alt)) {
      record.setAlt(alt);
    }
    if (isDefined(qual)) {
      record.setQual(Integer.parseInt(qual));
    }
    if (isDefined(filter)) {
      record.setFilter(filter);
    }
    if (isDefined(info)) {
      record.setInfo(info);
    }
    if (isDefined(format)) {
      record.setFormat(format);
    }
    for (int i = 9; i < fields.length; ++i) {
      record.addSample(fields[i]);
    }
    return record.build();
  }

  private static boolean isDefined(String value) {
    return !(null == value || ".".equals(value));
  }

  private final BufferedReader in;

  private VcfReader(BufferedReader in) {
    this.in = in;
  }

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
          Map<String, String> map = parseCompoundValue(value);
          MetaInformation.Info.Builder info = MetaInformation.Info.builder()
              .setId(map.get("ID"))
              .setNumber(MetaInformation.Number.create(map.get("Number")))
              .setType(MetaInformation.Info.Type.parse(map.get("Type")))
              .setDescription(map.get("Description"));
          for (Map.Entry<String, String> entry :
              getExtraFields(map, MetaInformation.Info.REQUIRED_FIELDS)) {
            info.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addInfo(info);
          break;
        case "FILTER":
          MetaInformation.Filter.Builder filter = MetaInformation.Filter.builder()
              .setId((map = parseCompoundValue(value)).get("ID"))
              .setDescription(map.get("Description"));
          for (Map.Entry<String, String> entry :
              getExtraFields(map, MetaInformation.Filter.REQUIRED_FIELDS)) {
            filter.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addFilter(filter);
          break;
        case "FORMAT":
          MetaInformation.Format.Builder format = MetaInformation.Format.builder()
              .setId((map = parseCompoundValue(value)).get("ID"))
              .setNumber(MetaInformation.Number.create(map.get("Number")))
              .setType(MetaInformation.Format.Type.parse(map.get("Type")))
              .setDescription(map.get("Description"));
          for (Map.Entry<String, String> entry :
              getExtraFields(map, MetaInformation.Format.REQUIRED_FIELDS)) {
            format.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addFormat(format);
          break;
        case "ALT":
          MetaInformation.Alt.Builder alt = MetaInformation.Alt.builder()
              .setId(MetaInformation.Alt.Type.parse((map = parseCompoundValue(value)).get("ID")))
              .setDescription(map.get("Description"));
          for (Map.Entry<String, String> entry :
              getExtraFields(map, MetaInformation.Alt.REQUIRED_FIELDS)) {
            alt.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addAlt(alt);
          break;
        case "assembly":
          metainfo.addAssembly(MetaInformation.Assembly.create(url(value)));
          break;
        case "contig":
          MetaInformation.Contig.Builder contig = MetaInformation.Contig.builder()
              .setId((map = parseCompoundValue(value, false)).get("ID"));
          for (Map.Entry<String, String> entry :
              getExtraFields(map, MetaInformation.Contig.REQUIRED_FIELDS)) {
            contig.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addContig(contig);
          break;
        case "SAMPLE":
          metainfo.addSample(MetaInformation.Sample.builder()
              .setId((map = parseCompoundValue(value)).get("ID"))
              .setGenome(map.get("Genome"))
              .setMixture(map.get("Mixture"))
              .setDescription(map.get("Description"))
              .build());
          break;
        case "PEDIGREE":
          MetaInformation.Pedigree.Builder pedigree = MetaInformation.Pedigree.builder();
          for (Map.Entry<String, String> entry : parseCompoundValue(value).entrySet()) {
            pedigree.addExtraField(entry.getKey(), entry.getValue());
          }
          metainfo.addPedigree(pedigree);
          break;
        case "pedigreeDB":
          metainfo.addPedigreeDB(MetaInformation.PedigreeDB.create(url(value)));
          break;
        default:
          metainfo.addUnparsedLine(MetaInformation.UnparsedMetaInfoLine.create(type, value));
      }
    }
    return metainfo.build();
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

  private static Map<String, String> parseCompoundValue(String value) {
    return parseCompoundValue(value, true);
  }

  private static Iterable<Map.Entry<String, String>> getExtraFields(Map<String, String> map, Set<String> requiredFields) {
    Map<String, String> extraFields = Maps.newLinkedHashMap(map);
    for (String field : MetaInformation.Info.REQUIRED_FIELDS) {
      extraFields.remove(field);
    }
    return extraFields.entrySet();
  }

  private static URL url(String url) {
    try {
      return new URL(url);
    } catch (MalformedURLException e) {
      throw Throwables.propagate(e);
    }
  }

  private static Map<String, String> parseCompoundValue(String input, boolean unescape) {
    if (input.startsWith("<") && input.endsWith(">")) {
      ImmutableMap.Builder<String, String> map = ImmutableMap.builder();
      for (
          Matcher matcher = KEY_VALUE_PAIR_PATTERN.matcher(input.substring(1, input.length() - 1));
          matcher.find(); ) {
        String value = matcher.group(2);
        map.put(matcher.group(1), unescape && value.startsWith("\"") && value.endsWith("\"")
            ? UNESCAPE_PATTERN.matcher(value.substring(1, value.length() - 1)).replaceAll("$1")
            : value);
        if (matcher.group(3).isEmpty()) {
          break;
        }
      }
      return map.build();
    }
    throw new IllegalStateException(String.format("Malformed value: \"%s\"", input));
  }
}

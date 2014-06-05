package edu.berkeley.cs.amplab.fastaparser;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Ordering;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.Comparator;
import java.util.Objects;
import java.util.SortedSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FastaIndex {

  public static final class Entry implements Comparable<Entry> {

    private static final Comparator<Entry> COMPARATOR = Ordering.natural()
        .onResultOf(
            new Function<Entry, Long>() {
              @Override public Long apply(Entry entry) {
                return entry.offset();
              }
            })
        .compound(Ordering.natural()
            .onResultOf(
                new Function<Entry, String>() {
                  @Override public String apply(Entry entry) {
                    return entry.name();
                  }
                }))
        .compound(Ordering.natural()
            .onResultOf(
                new Function<Entry, Long>() {
                  @Override public Long apply(Entry entry) {
                    return entry.length();
                  }
                }))
        .compound(Ordering.natural()
            .onResultOf(
                new Function<Entry, Integer>() {
                  @Override public Integer apply(Entry entry) {
                    return entry.bases();
                  }
                }))
        .compound(Ordering.natural()
            .onResultOf(
                new Function<Entry, Integer>() {
                  @Override public Integer apply(Entry entry) {
                    return entry.bytes();
                  }
                }));

    private static final Pattern PATTERN = Pattern.compile(
        "([^\t]+?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)");

    public static Entry parse(String line) {
      Matcher matcher = PATTERN.matcher(line);
      if (!matcher.matches()) {
        throw new IllegalArgumentException(
            String.format("Couldn't match \"%s\" against \"%s\"", line, PATTERN.pattern()));
      }
      return create(
          matcher.group(1),
          Long.parseLong(matcher.group(2)),
          Long.parseLong(matcher.group(3)),
          Integer.parseInt(matcher.group(4)),
          Integer.parseInt(matcher.group(5)));
    }

    public static Entry create(String name, long length, long offset, int bases, int bytes) {
      return new Entry(name, length, offset, bases, bytes);
    }

    private final String name;
    private final long length;
    private final long offset;
    private final int bases;
    private final int bytes;

    private Entry(String name, long length, long offset, int bases, int bytes) {
      this.name = name;
      this.length = length;
      this.offset = offset;
      this.bases = bases;
      this.bytes = bytes;
    }

    public int bases() {
      return bases;
    }

    public int bytes() {
      return bytes;
    }

    @Override
    public int compareTo(Entry rhs) {
      return COMPARATOR.compare(this, rhs);
    }

    @Override
    public boolean equals(Object obj) {
      return this == obj
          || null != obj
          && Entry.class == obj.getClass()
          && 0 == compareTo((Entry) obj);
    }

    @Override
    public int hashCode() {
      return Objects.hash(name(), length(), offset(), bases(), bytes());
    }

    public long length() {
      return length;
    }

    public String name() {
      return name;
    }

    public long offset() {
      return offset;
    }

    @Override
    public String toString() {
      return String.format("%s\t%d\t%d\t%d\t%d", name(), length(), offset(), bases(), bytes());
    }
  }

  public static FastaIndex create(File file) throws IOException {
    try (Reader in = new FileReader(file)) {
      return create(in);
    }
  }

  public static FastaIndex create(InputStream in) throws IOException {
    return create(new InputStreamReader(in));
  }

  public static FastaIndex create(Reader reader) throws IOException {
    @SuppressWarnings("resource") BufferedReader in = reader instanceof BufferedReader
        ? (BufferedReader) reader
        : new BufferedReader(reader);
    ImmutableSortedSet.Builder<Entry> entries = ImmutableSortedSet.naturalOrder();
    for (String line = in.readLine(); null != line; line = in.readLine()) {
      entries.add(Entry.parse(line));
    }
    return create(entries.build());
  }

  public static FastaIndex create(SortedSet<Entry> entries) {
    return new FastaIndex(entries);
  }

  private final SortedSet<Entry> entries;

  private FastaIndex(SortedSet<Entry> entries) {
    this.entries = entries;
  }

  @Override
  public boolean equals(Object obj) {
    return this == obj
        || null != obj
        && FastaIndex.class == obj.getClass()
        && Objects.equals(entries, ((FastaIndex) obj).entries);
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(entries);
  }

  @Override
  public String toString() {
    return String.format("%s%n", Joiner.on(String.format("%n")).join(entries));
  }

  public FastaIndex write(File file) throws IOException {
    try (Writer out = new FileWriter(file)) {
      return write(out);
    }
  }

  public FastaIndex write(OutputStream out) throws IOException {
    return write(new OutputStreamWriter(out));
  }

  public FastaIndex write(Writer writer) throws IOException {
    PrintWriter out = writer instanceof PrintWriter
        ? (PrintWriter) writer
        : new PrintWriter(writer);
    out.print(this);
    if (out.checkError()) {
      throw new IOException();
    }
    return this;
  }
}

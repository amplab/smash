/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.berkeley.cs.amplab.calldiff;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

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
import java.nio.CharBuffer;
import java.nio.charset.StandardCharsets;
import java.util.Comparator;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * An abstraction for parsing and making use of FASTA index files.
 */
public class FastaIndex {

  public static final class Entry implements Comparable<Entry> {

    private static final Comparator<Entry> COMPARATOR = Comparator
        .comparing(Entry::offset)
        .thenComparing(Entry::name)
        .thenComparing(Entry::length)
        .thenComparing(Entry::bases)
        .thenComparing(Entry::bytes);

    private static final Pattern PATTERN = Pattern.compile(
        "([^\t]+?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)\t(0|[1-9][0-9]*?)");

    public static Entry create(String name, long length, long offset, int bases, int bytes) {
      return new Entry(name, length, offset, bases, bytes);
    }

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

    private final int bases;
    private final int bytes;
    private final long length;
    private final String name;
    private final long offset;

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
      return Objects.hash(fields().toArray());
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
      return fields().map(Object::toString).collect(Collectors.joining("\t"));
    }

    private Stream<Object> fields() {
      return Stream.of(name(), length(), offset(), bases(), bytes());
    }
  }

  private static final int NEWLINE_SIZE_IN_BYTES =
      StandardCharsets.UTF_8.encode(CharBuffer.wrap(String.format("%n"))).limit();

  public static FastaIndex create(File fastaFile) throws IOException {
    class ExceptionWrapper extends RuntimeException {

      ExceptionWrapper(IOException cause) {
        super(cause);
      }

      @Override public IOException getCause() {
        return (IOException) super.getCause();
      }
    }
    try (final BufferedReader in = new BufferedReader(new FileReader(fastaFile))) {
      ImmutableSortedSet.Builder<FastaIndex.Entry> entries = ImmutableSortedSet.naturalOrder();
      String name = null;
      long length = 0;
      long offset1 = 0;
      long offset2 = 0;
      int bases = 0;
      for (
          PeekingIterator<String> iterator = Iterators.peekingIterator(
              new AbstractIterator<String>() {
                @Override protected String computeNext() {
                  try {
                    String line = in.readLine();
                    return null == line ? endOfData() : line;
                  } catch (IOException e) {
                    throw new ExceptionWrapper(e);
                  }
                }
              });
          iterator.hasNext();) {
        String line = iterator.next();
        int lineLength = line.length();
        offset2 += lineLength + NEWLINE_SIZE_IN_BYTES;
        if (line.startsWith(">")) {
          if (null != name) {
            entries.add(FastaIndex.Entry.create(
                name, length, offset1, bases, bases + NEWLINE_SIZE_IN_BYTES));
          }
          name = line.substring(1).split("\\p{Space}+?")[0];
          length = 0;
          offset1 = offset2;
          bases = 0;
        } else {
          if (0 == bases) {
            bases = lineLength;
          } else if (bases != lineLength
              && iterator.hasNext()
              && !iterator.peek().startsWith(">")) {
            throw new IllegalStateException(String.format(
                "Inconsistent line lengths at contig \"%s\", offset %d",
                name,
                offset2 - lineLength - NEWLINE_SIZE_IN_BYTES));
          }
          length += lineLength;
        }
      }
      if (null != name) {
        entries.add(FastaIndex.Entry.create(
            name, length, offset1, bases, bases + NEWLINE_SIZE_IN_BYTES));
      }
      return createFromEntries(entries.build());
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }

  public static FastaIndex createFromEntries(SortedSet<Entry> entries) {
    return new FastaIndex(entries);
  }

  public static FastaIndex read(File file) throws IOException {
    try (Reader in = new FileReader(file)) {
      return read(in);
    }
  }

  public static FastaIndex read(InputStream in) {
    return read(new InputStreamReader(in));
  }

  public static FastaIndex read(Reader reader) {
    @SuppressWarnings("resource") BufferedReader in = reader instanceof BufferedReader
        ? (BufferedReader) reader
        : new BufferedReader(reader);
    return createFromEntries(
        in.lines().map(Entry::parse).collect(Collectors.toCollection(TreeSet::new)));
  }

  private final SortedSet<Entry> entries;

  private FastaIndex(SortedSet<Entry> entries) {
    this.entries = entries;
  }

  public SortedSet<Entry> entries() {
    return entries;
  }

  @Override
  public boolean equals(Object obj) {
    return this == obj
        || null != obj
        && FastaIndex.class == obj.getClass()
        && Objects.equals(entries(), ((FastaIndex) obj).entries());
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(entries());
  }

  @Override
  public String toString() {
    return entries().stream()
        .map(Object::toString)
        .collect(Collectors.joining(String.format("%n")));
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
    @SuppressWarnings("resource") PrintWriter out = writer instanceof PrintWriter
        ? (PrintWriter) writer
        : new PrintWriter(writer);
    out.println(this);
    if (out.checkError()) {
      throw new IOException();
    }
    return this;
  }
}

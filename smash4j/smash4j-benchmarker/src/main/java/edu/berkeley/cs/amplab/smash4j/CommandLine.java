package edu.berkeley.cs.amplab.smash4j;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;
import com.google.common.base.Optional;

import java.io.File;
import java.io.PrintStream;

@Parameters(separators = "=")
public class CommandLine {

  public static Optional<CommandLine> parse(String... args) {
    CommandLine commandLine = new CommandLine();
    JCommander jCommander = new JCommander(commandLine);
    jCommander.setProgramName("smash4j");
    try {
      jCommander.parse(args);
      if (commandLine.help) {
        return usage(System.out, jCommander);
      }
      return Optional.of(commandLine);
    } catch (ParameterException e) {
      return usage(System.err.format("%s%n%n", e.getMessage()), jCommander);
    }
  }

  private static Optional<CommandLine> usage(PrintStream out, JCommander jCommander) {
    StringBuilder buffer = new StringBuilder();
    jCommander.usage(buffer);
    out.print(buffer);
    return Optional.absent();
  }

  @Parameter(
      names = { "-i", "--in" },
      description = "Input file",
      required = true,
      converter = FileConverter.class)
  private File in;

  @Parameter(
      names = { "-o", "--out" },
      description = "Output file",
      required = true,
      converter = FileConverter.class)
  private File out;

  @Parameter(
      names = { "-h", "--help" },
      description = "Print this help message",
      help = true)
  private boolean help;

  public File in() {
    return in;
  }

  public File out() {
    return out;
  }
}

package edu.berkeley.cs.amplab.vardiff;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.Optional;

@Parameters(separators = "=")
public class CommandLine {

  public static class Builder {

    private String lhsSampleId;
    private String lhsVcf;
    private String referenceFai;
    private String referenceFasta;
    private String rhsSampleId;
    private String rhsVcf;
    private boolean presorted;

    public CommandLine build() {
      return new CommandLine(
          lhsSampleId,
          lhsVcf,
          referenceFai,
          referenceFasta,
          rhsSampleId,
          rhsVcf,
          presorted);
    }

    public Builder setLhsSampleId(String lhsSampleId) {
      this.lhsSampleId = lhsSampleId;
      return this;
    }

    public Builder setLhsVcf(String lhsVcf) {
      this.lhsVcf = lhsVcf;
      return this;
    }

    public Builder setPresorted(boolean presorted) {
      this.presorted = presorted;
      return this;
    }

    public Builder setReferenceFai(String referenceFai) {
      this.referenceFai = referenceFai;
      return this;
    }

    public Builder setReferenceFasta(String referenceFasta) {
      this.referenceFasta = referenceFasta;
      return this;
    }

    public Builder setRhsSampleId(String rhsSampleId) {
      this.rhsSampleId = rhsSampleId;
      return this;
    }

    public Builder setRhsVcf(String rhsVcf) {
      this.rhsVcf = rhsVcf;
      return this;
    }
  }

  private static final HashCodeAndEquals<CommandLine>
      HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
          CommandLine.class,
          CommandLine::lhsSampleId,
          CommandLine::lhsVcf,
          CommandLine::referenceFai,
          CommandLine::referenceFasta,
          CommandLine::rhsSampleId,
          CommandLine::rhsVcf,
          CommandLine::presorted);

  public static Builder builder() {
    return new Builder();
  }

  public static CommandLine parse(String... args) {
    CommandLine commandLine = new CommandLine();
    new JCommander(commandLine, args);
    return commandLine;
  }

  @Parameter(names = { "--lhs_sample_id" })
  private String lhsSampleId;

  @Parameter(names = { "--lhs_vcf" })
  private String lhsVcf;

  @Parameter(names = { "--reference_fai" })
  private String referenceFai;

  @Parameter(names = { "--reference_fasta" })
  private String referenceFasta;

  @Parameter(names = { "--rhs_sample_id" })
  private String rhsSampleId;

  @Parameter(names = { "--rhs_vcf" })
  private String rhsVcf;

  @Parameter(names = { "--presorted" })
  private boolean presorted;

  private CommandLine() {
    this(null, null, null, null, null, null, false);
  }

  private CommandLine(
      String lhsSampleId,
      String lhsVcf,
      String referenceFai,
      String referenceFasta,
      String rhsSampleId,
      String rhsVcf,
      boolean presorted) {
    this.lhsSampleId = lhsSampleId;
    this.lhsVcf = lhsVcf;
    this.referenceFai = referenceFai;
    this.referenceFasta = referenceFasta;
    this.rhsSampleId = rhsSampleId;
    this.rhsVcf = rhsVcf;
    this.presorted = presorted;
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public Optional<String> lhsSampleId() {
    return Optional.ofNullable(lhsSampleId);
  }

  public Optional<String> lhsVcf() {
    return Optional.ofNullable(lhsVcf);
  }

  public boolean presorted() {
    return presorted;
  }

  public Optional<String> referenceFai() {
    return Optional.ofNullable(referenceFai);
  }

  public Optional<String> referenceFasta() {
    return Optional.ofNullable(referenceFasta);
  }

  public Optional<String> rhsSampleId() {
    return Optional.ofNullable(rhsSampleId);
  }

  public Optional<String> rhsVcf() {
    return Optional.ofNullable(rhsVcf);
  }
}

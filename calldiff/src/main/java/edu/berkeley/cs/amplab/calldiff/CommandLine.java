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

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.Optional;

/**
 * The class responsible for parsing and storing the command line used to invoke the program.
 */
@Parameters(separators = "=")
public class CommandLine {

  public static class Builder {

    private String apiKey;
    private String clientSecretsFile;
    private String lhsCallsetId;
    private String lhsSampleId;
    private String lhsVcf;
    private String p12File;
    private boolean presorted;
    private String referenceFai;
    private String referenceFasta;
    private String rhsCallsetId;
    private String rhsSampleId;
    private String rhsVcf;
    private String rootUrl;
    private String serviceAccountId;
    private Integer timeout;

    public CommandLine build() {
      return new CommandLine(
          apiKey,
          clientSecretsFile,
          lhsCallsetId,
          lhsSampleId,
          lhsVcf,
          p12File,
          presorted,
          referenceFai,
          referenceFasta,
          rhsCallsetId,
          rhsSampleId,
          rhsVcf,
          rootUrl,
          serviceAccountId,
          timeout);
    }

    public Builder setApiKey(String apiKey) {
      this.apiKey = apiKey;
      return this;
    }

    public Builder setClientSecretsFile(String clientSecretsFile) {
      this.clientSecretsFile = clientSecretsFile;
      return this;
    }

    public Builder setLhsCallsetId(String lhsCallsetId) {
      this.lhsCallsetId = lhsCallsetId;
      return this;
    }

    public Builder setLhsSampleId(String lhsSampleId) {
      this.lhsSampleId = lhsSampleId;
      return this;
    }

    public Builder setLhsVcf(String lhsVcf) {
      this.lhsVcf = lhsVcf;
      return this;
    }

    public Builder setP12File(String p12File) {
      this.p12File = p12File;
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

    public Builder setRhsCallsetId(String rhsCallsetId) {
      this.rhsCallsetId = rhsCallsetId;
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

    public Builder setRootUrl(String rootUrl) {
      this.rootUrl = rootUrl;
      return this;
    }

    public Builder setServiceAccountId(String serviceAccountId) {
      this.serviceAccountId = serviceAccountId;
      return this;
    }

    public Builder setTimeout(Integer timeout) {
      this.timeout = timeout;
      return this;
    }
  }

  private static final HashCodeAndEquals<CommandLine>
      HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
          CommandLine.class,
          CommandLine::apiKey,
          CommandLine::clientSecretsFile,
          CommandLine::lhsCallsetId,
          CommandLine::lhsSampleId,
          CommandLine::lhsVcf,
          CommandLine::p12File,
          CommandLine::presorted,
          CommandLine::referenceFai,
          CommandLine::referenceFasta,
          CommandLine::rhsCallsetId,
          CommandLine::rhsSampleId,
          CommandLine::rhsVcf,
          CommandLine::rootUrl,
          CommandLine::serviceAccountId,
          CommandLine::timeout);

  public static Builder builder() {
    return new Builder();
  }

  public static Optional<CommandLine> parse(String... args) {
    CommandLine commandLine = new CommandLine();
    JCommander jCommander = new JCommander(commandLine);
    jCommander.setProgramName("java -jar calldiff-jar-with-dependencies.jar");
    jCommander.parse(args);
    if (commandLine.help) {
      StringBuilder buffer = new StringBuilder();
      jCommander.usage(buffer);
      System.err.print(buffer);
      return Optional.empty();
    }
    return Optional.of(commandLine);
  }

  @Parameter(
      names = { "--api_key" },
      description = "The API key used to authenticate to your Google Cloud project")
  private String apiKey;

  @Parameter(
      names = { "--client_secrets_file" },
      description = "The client secrets file used to authorize access to your Google Cloud project")
  private String clientSecretsFile;

  @Parameter(
      names = { "--help" },
      description = "Print the help message",
      help = true)
  private boolean help;

  @Parameter(
      names = { "--lhs_callset_id" },
      description = "The callset id to use on the left hand side of the comparison")
  private String lhsCallsetId;

  @Parameter(
      names = { "--lhs_sample_id" },
      description = "The sample id to use on the left hand side of the comparison")
  private String lhsSampleId;

  @Parameter(
      names = { "--lhs_vcf" },
      description = "The path to the VCF file to use on the left hand side of the comparison")
  private String lhsVcf;

  @Parameter(
      names = { "--p12_file" },
      description = "The P12 file containing the private key that authorizes the service account "
          + "for your Google Cloud Project")
  private String p12File;

  @Parameter(
      names = { "--presorted" },
      description = "Skip sorting the input because it is already properly sorted")
  private boolean presorted;

  @Parameter(
      names = { "--reference_fai" },
      description = "The FASTA index file for the reference sequence")
  private String referenceFai;

  @Parameter(
      names = { "--reference_fasta" },
      description = " The FASTA file for the reference sequence")
  private String referenceFasta;

  @Parameter(
      names = { "--rhs_callset_id" },
      description = "The callset id to use on the right hand side of the comparison")
  private String rhsCallsetId;

  @Parameter(
      names = { "--rhs_sample_id" },
      description = "The sample id to use on the right hand side of the comparison")
  private String rhsSampleId;

  @Parameter(
      names = { "--rhs_vcf" },
      description = "The path to the VCF file to use on the right hand side of the comparison")
  private String rhsVcf;

  @Parameter(
      names = { "--root_url" },
      description = "The URL to communicate with to fetch variants from the cloud")
  private String rootUrl;

  @Parameter(
      names = { "--service_account_id" },
      description = "The email address for the service account used to authorize your Google "
          + "Cloud project")
  private String serviceAccountId;

  @Parameter(
      names = { "--timeout" },
      description = "The connect and read timeouts to use when making requests to the cloud")
  private Integer timeout;

  public CommandLine() {
    this(null, null, null, null, null, null, false, null, null, null, null, null, null, null, null);
  }

  private CommandLine(
      String apiKey,
      String clientSecretsFile,
      String lhsCallsetId,
      String lhsSampleId,
      String lhsVcf,
      String p12File,
      boolean presorted,
      String referenceFai,
      String referenceFasta,
      String rhsCallsetId,
      String rhsSampleId,
      String rhsVcf,
      String rootUrl,
      String serviceAccountId,
      Integer timeout) {
    this.apiKey = apiKey;
    this.clientSecretsFile = clientSecretsFile;
    this.lhsCallsetId = lhsCallsetId;
    this.lhsSampleId = lhsSampleId;
    this.lhsVcf = lhsVcf;
    this.p12File = p12File;
    this.presorted = presorted;
    this.referenceFai = referenceFai;
    this.referenceFasta = referenceFasta;
    this.rhsCallsetId = rhsCallsetId;
    this.rhsSampleId = rhsSampleId;
    this.rhsVcf = rhsVcf;
    this.rootUrl = rootUrl;
    this.serviceAccountId = serviceAccountId;
    this.timeout = timeout;
  }

  public Optional<String> apiKey() {
    return Optional.ofNullable(apiKey);
  }

  public Optional<String> clientSecretsFile() {
    return Optional.ofNullable(clientSecretsFile);
  }

  @Override
  public boolean equals(Object obj) {
    return HASH_CODE_AND_EQUALS.equals(this, obj);
  }

  @Override
  public int hashCode() {
    return HASH_CODE_AND_EQUALS.hashCode(this);
  }

  public Optional<String> lhsCallsetId() {
    return Optional.ofNullable(lhsCallsetId);
  }

  public Optional<String> lhsSampleId() {
    return Optional.ofNullable(lhsSampleId);
  }

  public Optional<String> lhsVcf() {
    return Optional.ofNullable(lhsVcf);
  }

  public Optional<String> p12File() {
    return Optional.ofNullable(p12File);
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

  public Optional<String> rhsCallsetId() {
    return Optional.ofNullable(rhsCallsetId);
  }

  public Optional<String> rhsSampleId() {
    return Optional.ofNullable(rhsSampleId);
  }

  public Optional<String> rhsVcf() {
    return Optional.ofNullable(rhsVcf);
  }

  public Optional<String> rootUrl() {
    return Optional.ofNullable(rootUrl);
  }

  public Optional<String> serviceAccountId() {
    return Optional.ofNullable(serviceAccountId);
  }

  public Optional<Integer> timeout() {
    return Optional.ofNullable(timeout);
  }
}

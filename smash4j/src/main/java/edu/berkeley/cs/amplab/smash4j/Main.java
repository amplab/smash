package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.smash4j.CommandDispatcher.MainCommand;
import edu.berkeley.cs.amplab.smash4j.fasta.FastaReader;

import java.io.File;
import java.util.Map;
import java.util.prefs.Preferences;

public class Main {

  static final Preferences
      PREFERENCES = Preferences.userNodeForPackage(Main.class);
  static final String
      PREFERENCES_PATH = PREFERENCES.absolutePath(),
      PROGRAM_NAME = "AMPLab-SMaSH4J/0.1";

  private static void main(final MainCommand command) throws Exception {
    File reference = command.referenceFasta();
    Optional<File> index = command.referenceFastaIndex();
    output(
        (index.isPresent() ? FastaReader.create(reference, index.get()) : FastaReader.create(reference)).read(
            new FastaReader.Callback<Map<String, VariantEvaluator.ContigStats>>() {

              private VariantScanner createScanner(Optional<String> callset, Optional<File> vcf) {
                if (callset.isPresent() && !vcf.isPresent()) {
                  return VariantScanner.fromCallset(callset.get());
                } else if (!callset.isPresent() && vcf.isPresent()) {
                  return VariantScanner.fromVcfFile(vcf.get());
                } else {
                  throw new IllegalStateException();
                }
              }

              @Override public Map<String, VariantEvaluator.ContigStats> read(
                  Map<String, Integer> info,
                  final FastaReader.FastaFile reference) throws Exception {
                return createScanner(command.trueCallset(), command.trueVcf()).scan(
                    new VariantScanner.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                      @Override public Map<String, VariantEvaluator.ContigStats> scan(
                          final FluentIterable<Variant> trueVariants) throws Exception {
                        return createScanner(command.predictedCallset(), command.predictedVcf()).scan(
                            new VariantScanner.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                              @Override public Map<String, VariantEvaluator.ContigStats> scan(
                                  FluentIterable<Variant> predictedVariants) throws Exception {
                                int maxIndelSize = command.maxIndelSize().or(50);
                                Normalizer normalizer = Normalizer.create(maxIndelSize, reference);
                                return VariantEvaluator.builder()
                                    .setMaxIndexSize(maxIndelSize)
                                    .setMaxSvBreakpointDistance(
                                        command.maxSvBreakpointDistance().or(100))
                                    .setMaxVariantLengthDifference(
                                        command.maxVariantLengthDifference().or(100))
                                    .setReference(reference)
                                    .setRescueWindowSize(command.rescueWindowSize().or(50))
                                    .build()
                                    .evaluate(
                                        normalizer.normalize(trueVariants),
                                        normalizer.normalize(predictedVariants));
                              }
                            });
                      }
                    });
              }
            }));
  }

  private static void output(Map<String, VariantEvaluator.ContigStats> contigStats) {
    // TODO
  }

  public static void main(String[] args) throws Exception {
    new CommandDispatcher() {

      @Override protected void main(MainCommand command) throws Exception {
        Main.main(command);
      }

      @Override protected void setPrefs(CommandDispatcher.SetPrefsCommand command) {
        Optional<CommandDispatcher.SetPrefsCommand.AuthorizationMethod>
            authorizationMethod = command.authorizationMethod();
        if (authorizationMethod.isPresent()) {
          PREFERENCES.put("authorizationMethod", authorizationMethod.get().toString());
        }
        Optional<String> apiKey = command.apiKey();
        if (apiKey.isPresent()) {
          PREFERENCES.put("apiKey", apiKey.get());
        }
        Optional<File> clientSecretsFile = command.clientSecretsFile();
        if (clientSecretsFile.isPresent()) {
          PREFERENCES.put("clientSecretsFile", clientSecretsFile.get().getPath());
        }
        Optional<String> serviceAccountId = command.serviceAccountId();
        if (serviceAccountId.isPresent()) {
          PREFERENCES.put("serviceAccountId", serviceAccountId.get());
        }
        Optional<File> serviceAccountP12File = command.serviceAccountP12File();
        if (serviceAccountP12File.isPresent()) {
          PREFERENCES.put("serviceAccountP12File", serviceAccountP12File.get().getPath());
        }
      }

      private void showPref(boolean condition, String key) {
        if (condition) {
          String value = PREFERENCES.get(key, null);
          System.out.format("%s:%s = %s%n",
              PREFERENCES_PATH, key, null == value ? null : String.format("\"%s\"", value));
        }
      }

      @Override protected void showPrefs(CommandDispatcher.ShowPrefsCommand command) {
        boolean
            authorizationMethod = command.authorizationMethod(),
            apiKey = command.apiKey(),
            clientSecretsFile = command.clientSecretsFile(),
            serviceAccountId = command.serviceAccountId(),
            serviceAccountP12File = command.serviceAccountP12File(),
            noFlags = !(authorizationMethod
                || apiKey
                || clientSecretsFile
                || serviceAccountId
                || serviceAccountP12File);
        showPref(authorizationMethod || noFlags, "authorizationMethod");
        showPref(apiKey || noFlags, "apiKey");
        showPref(clientSecretsFile || noFlags, "clientSecretsFile");
        showPref(serviceAccountId || noFlags, "serviceAccountId");
        showPref(serviceAccountP12File || noFlags, "serviceAccountP12File");
      }
    }.parse(args);
  }
}

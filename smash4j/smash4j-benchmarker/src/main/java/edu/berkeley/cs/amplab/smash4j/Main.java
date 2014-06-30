package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Multimap;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;
import edu.berkeley.cs.amplab.smash4j.VariantEvaluator.ContigStats;
import edu.berkeley.cs.amplab.smash4j.VariantEvaluator.GenotypeConcordance;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Map;
import java.util.NavigableMap;
import java.util.prefs.Preferences;

public class Main {

  static final Preferences
      PREFERENCES = Preferences.userNodeForPackage(Main.class);
  static final String
      PREFERENCES_PATH = PREFERENCES.absolutePath(),
      PROGRAM_NAME = "AMPLab-SMaSH4J/0.1";

  public static void main(String[] args) throws Exception {
    new CommandDispatcher() {

      private Optional<VariantScanner> createScanner(Optional<String> callset, Optional<File> vcf) {
        return callset.isPresent() && !vcf.isPresent()
            ? Optional.of(VariantScanner.fromCallsets(Collections.singletonList(callset.get())))
            : !callset.isPresent() && vcf.isPresent()
                ? Optional.of(VariantScanner.fromVcfFile(vcf.get()))
                : Optional.<VariantScanner>absent();
      }

      @Override protected void diff(final DiffCommand command) throws Exception {
        File fasta = command.referenceFasta();
        Optional<File> index = command.referenceFastaIndex();
        print(System.out, (index.isPresent() ? FastaReader.create(fasta, index.get()) : FastaReader.create(fasta))
            .read(
                new FastaReader.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                  @Override public Map<String, VariantEvaluator.ContigStats> read(Map<String, Integer> info,
                      final FastaReader.Callback.FastaFile reference) throws Exception {
                    return createScanner(command.trueCallset(), command.trueVcf()).get().scan(
                        new VariantScanner.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                          @Override public Map<String, VariantEvaluator.ContigStats> scan(
                              final FluentIterable<VariantProto> trueVars) throws Exception {
                            return createScanner(command.predictedCallset(), command.predictedVcf()).get().scan(
                                new VariantScanner.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                                  @Override public Map<String, VariantEvaluator.ContigStats> scan(
                                      final FluentIterable<VariantProto> predictedVars) throws Exception {
                                    VariantEvaluator.Builder builder = VariantEvaluator.builder()
                                        .setReference(reference);
                                    Optional<Integer>
                                        maxIndelSize = command.maxIndelSize(),
                                        maxSvBreakpointDist = command.maxSvBreakpointDistance(),
                                        maxVariantLengthDiff = command.maxVariantLengthDifference(),
                                        rescueWindowSize = command.rescueWindowSize();
                                    if (maxIndelSize.isPresent()) {
                                      builder.setMaxIndexSize(
                                          maxIndelSize.get());
                                    }
                                    if (maxSvBreakpointDist.isPresent()) {
                                      builder.setMaxSvBreakpointDistance(
                                          maxSvBreakpointDist.get());
                                    }
                                    if (maxVariantLengthDiff.isPresent()) {
                                      builder.setMaxVariantLengthDifference(
                                          maxVariantLengthDiff.get());
                                    }
                                    if (rescueWindowSize.isPresent()) {
                                      builder.setRescueWindowSize(
                                          rescueWindowSize.get());
                                    }
                                    final VariantEvaluator evaluator = builder.build();
                                    Optional<VariantScanner> fpScanner = createScanner(
                                        command.knownFalsePositivesCallset(),
                                        command.knownFalsePositivesVcf());
                                    return fpScanner.isPresent()
                                        ? fpScanner.get().scan(
                                            new VariantScanner.Callback<Map<String, VariantEvaluator.ContigStats>>() {
                                              @Override public Map<String, ContigStats> scan(
                                                  FluentIterable<VariantProto> knownFalsePositives) {
                                                return evaluator.evaluate(
                                                    trueVars, predictedVars, knownFalsePositives);
                                              }
                                            })
                                        : evaluator.evaluate(trueVars, predictedVars);
                                  }
                                });
                          }
                        });
                  }
                }));
      }

      @Override protected void main(MainCommand command) throws Exception {}

      @Override
      protected void normalize(final NormalizeCommand command) throws Exception {
        File fasta = command.referenceFasta();
        Optional<File> index = command.referenceFastaIndex();
        (index.isPresent() ? FastaReader.create(fasta, index.get()) : FastaReader.create(fasta))
            .read(new FastaReader.Callback<Void>() {
              @Override public Void read(Map<String, Integer> info,
                  final FastaReader.Callback.FastaFile fastaFile) throws Exception {
                return createScanner(command.callset(), command.vcf()).get().scan(
                    new VariantScanner.Callback<Void>() {
                      @Override
                      public Void scan(FluentIterable<VariantProto> variants) throws Exception {
                        try (OutputStream out = new FileOutputStream(command.out().get())) {
                          int maxIndelSize = command.maxIndexSize().or(50);
                          Normalizer normalizer = command.cleanOnly()
                              ? Normalizer.cleanOnly(maxIndelSize, fastaFile)
                              : Normalizer.create(maxIndelSize, fastaFile);
                          VariantWriter.create(out).writeVariants(
                              VariantScanner.fromVariantProtos(normalizer.normalize(variants)));
                          return null;
                        }
                      }
                    });
              }
            });
      }

      private PrintStream print(PrintStream out, Map<String, VariantEvaluator.ContigStats> contigStats) {
        for (VariantEvaluator.ContigStats stats : contigStats.values()) {
          Optional<NavigableMap<Integer, VariantProto>>
              allKnownFalsePositives = stats.allKnownFalsePositives(),
              correctKnownFalsePositives = stats.correctKnownFalsePositives(),
              knownFalsePositives = stats.knownFalsePositives(),
              rescuedVariants = stats.rescuedVariants();
          GenotypeConcordance
              concordance = stats.concordance();
          NavigableMap<Integer, VariantProto>
              falseNegatives = stats.falseNegatives(),
              falsePositives = stats.falsePositives(),
              predictedVariants = stats.predictedVariants(),
              truePositives = stats.truePositives(),
              trueVariants = stats.trueVariants();
          Multimap<VariantType, Integer>
              incorrectPredictions = stats.incorrectPredictions();
          
        }
        return out;
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

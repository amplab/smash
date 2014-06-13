package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;

import edu.berkeley.cs.amplab.fastaparser.FastaReader;
import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Collections;
import java.util.Map;
import java.util.prefs.Preferences;

public class Main {

  static final Preferences
      PREFERENCES = Preferences.userNodeForPackage(Main.class);
  static final String
      PREFERENCES_PATH = PREFERENCES.absolutePath(),
      PROGRAM_NAME = "AMPLab-SMaSH4J/0.1";

  public static void main(String[] args) throws Exception {
    new CommandDispatcher() {

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

      private void showPref(boolean condition, String key) {
        if (condition) {
          String value = PREFERENCES.get(key, null);
          System.out.format("%s:%s = %s%n",
              PREFERENCES_PATH, key, null == value ? null : String.format("\"%s\"", value));
        }
      }

      @Override
      protected void normalize(final NormalizeCommand command) throws Exception {
        File fasta = command.referenceFasta();
        Optional<File> index = command.referenceFastaIndex();
        (index.isPresent() ? FastaReader.create(fasta, index.get()) : FastaReader.create(fasta))
            .read(new FastaReader.Callback<Void>() {
              @Override public Void read(Map<String, Integer> info,
                  final FastaReader.Callback.FastaFile fastaFile) throws Exception {
                Optional<String> callset = command.callset();
                Optional<File> vcf = command.vcf();
                final VariantScanner scanner;
                if (callset.isPresent() && !vcf.isPresent()) {
                  scanner = VariantScanner.fromCallsets(Collections.singletonList(callset.get()));
                } else if (!callset.isPresent() && vcf.isPresent()) {
                  scanner = VariantScanner.fromVcfFile(vcf.get());
                } else {
                  throw new IllegalStateException();
                }
                return scanner.scan(
                    new VariantScanner.Callback<Void>() {
                      @Override
                      public Void scan(FluentIterable<VariantProto> variants) throws Exception {
                        try (OutputStream out = new FileOutputStream(command.out())) {
                          int maxIndelSize = command.maxIndexSize().or(50);
                          VariantWriter.create(out).writeVariants(VariantScanner.fromVariantProtos(
                              (command.cleanOnly()
                                  ? Normalizer.cleanOnly(maxIndelSize, fastaFile)
                                  : Normalizer.create(maxIndelSize, fastaFile))
                                  .normalize(variants)));
                          return null;
                        }
                      }
                    });
              }
            });
      }

      @Override protected void main(MainCommand command) throws Exception {
      }
    }.parse(args);
  }
}

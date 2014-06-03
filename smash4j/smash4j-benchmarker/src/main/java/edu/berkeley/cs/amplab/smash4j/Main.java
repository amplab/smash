package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableMap;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

import edu.berkeley.cs.amplab.smash4j.Smash4J.VariantProto;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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

      @Override protected void main(MainCommand command) throws Exception {
        Optional<String>
            leftHandSide = command.leftHandSide(),
            rightHandSide = command.rightHandSide();
        for (Map.Entry<String, Optional<String>> entry : ImmutableMap
            .of("lhs", leftHandSide, "rhs", rightHandSide)
            .entrySet()) {
          if (!entry.getValue().isPresent()) {
            throw new IllegalArgumentException(
                String.format("Flag --%s is required", entry.getKey()));
          }
        }
        Main.main(leftHandSide.get(), rightHandSide.get());
      }
    }.parse(args);
  }

  private static final ListeningExecutorService
      EXECUTOR_SERVICE = MoreExecutors.listeningDecorator(Executors.newCachedThreadPool());

  private static void main(String leftHandSide, String rightHandSide) throws Exception {
    File desktop = new File("/usr/local/google/home/kwestbrooks/Desktop");
    List<File> list = Futures.allAsList(
        cleanAndNormalize(leftHandSide, "lhs", desktop, false),
        cleanAndNormalize(rightHandSide, "rhs", desktop, false)).get();
    File lhs = list.get(0), rhs = list.get(1);
    
  }

  private static ListenableFuture<File> cleanAndNormalize(
      final String spec,
      final String prefix,
      final File directory,
      final boolean deleteOnExit) {
    return EXECUTOR_SERVICE.submit(
        new Callable<File>() {
          @Override public File call() throws Exception {
            return VariantScanner.create(spec).scan(
                new VariantScanner.Callback<File>() {
                  @Override public File accept(FluentIterable<VariantProto> variants)
                      throws Exception {
                    File tempFile = File.createTempFile(prefix, ".variants", directory);
                    if (deleteOnExit) {
                      tempFile.deleteOnExit();
                    }
                    try (OutputStream out = new FileOutputStream(tempFile)) {
                      for (VariantProto variant : variants) {
                        if (filter(variant)) {
                          normalize(variant).writeDelimitedTo(out);
                        }
                      }
                      return tempFile;
                    }
                  }
                });
          }
        });
  }

  private static boolean filter(VariantProto variant) {
    return true;
  }

  private static VariantProto normalize(VariantProto variant) {
    return variant;
  }
}

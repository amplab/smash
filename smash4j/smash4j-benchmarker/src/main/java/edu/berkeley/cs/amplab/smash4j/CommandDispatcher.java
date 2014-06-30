package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.IValueValidator;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.MissingCommandException;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

import java.io.File;

public abstract class CommandDispatcher {

  @Parameters(separators = "=")
  public static class DiffCommand {

    @Parameter(
        names = {"--indel_err"})
    private Double indelErr;

    @Parameter(
        names = { "--known_false_positives_callset" })
    private String knownFalsePositivesCallset;

    @Parameter(
        names = { "--known_false_positives_vcf" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File knownFalsePositivesVcf;

    @Parameter(
        names = { "--max_indel_size" })
    private Integer maxIndelSize;

    @Parameter(
        names = { "--max_sv_breakpoint_distance" })
    private Integer maxSvBreakpointDistance;

    @Parameter(
        names = { "--max_variant_length_difference" })
    private Integer maxVariantLengthDifference;

    @Parameter(
        names = { "--predicted_callset" })
    private String predictedCallset;

    @Parameter(
        names = { "--predicted_vcf" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File predictedVcf;

    @Parameter(
        names = { "--reference_fasta" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File referenceFasta;

    @Parameter(
        names = { "--reference_fasta_index" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File referenceFastaIndex;

    @Parameter(
        names = { "--rescue_window_size" })
    private Integer rescueWindowSize;

    @Parameter(
        names = {"--snp_err"})
    private Double snpErr;

    @Parameter(
        names = {"--sv_err"})
    private Double svErr;

    @Parameter(
        names = { "--true_callset" })
    private String trueCallset;

    @Parameter(
        names = { "--true_vcf" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File trueVcf;

    public Optional<Double> indelErr() {
      return Optional.fromNullable(indelErr);
    }

    public Optional<String> knownFalsePositivesCallset() {
      return Optional.fromNullable(knownFalsePositivesCallset);
    }

    public Optional<File> knownFalsePositivesVcf() {
      return Optional.fromNullable(knownFalsePositivesVcf);
    }

    public Optional<Integer> maxIndelSize() {
      return Optional.fromNullable(maxIndelSize);
    }

    public Optional<Integer> maxSvBreakpointDistance() {
      return Optional.fromNullable(maxSvBreakpointDistance);
    }

    public Optional<Integer> maxVariantLengthDifference() {
      return Optional.fromNullable(maxVariantLengthDifference);
    }

    public Optional<String> predictedCallset() {
      return Optional.fromNullable(predictedCallset);
    }

    public Optional<File> predictedVcf() {
      return Optional.fromNullable(predictedVcf);
    }

    public File referenceFasta() {
      return referenceFasta;
    }

    public Optional<File> referenceFastaIndex() {
      return Optional.fromNullable(referenceFastaIndex);
    }

    public Optional<Integer> rescueWindowSize() {
      return Optional.fromNullable(rescueWindowSize);
    }

    public Optional<Double> snpErr() {
      return Optional.fromNullable(snpErr);
    }

    public Optional<Double> svErr() {
      return Optional.fromNullable(svErr);
    }

    public Optional<String> trueCallset() {
      return Optional.fromNullable(trueCallset);
    }

    public Optional<File> trueVcf() {
      return Optional.fromNullable(trueVcf);
    }

    DiffCommand validate() {
      validateDouble(indelErr(), "indel_err");
      validateDouble(snpErr(), "snp_err");
      validateDouble(svErr(), "sv_err");
      validateInteger(maxIndelSize(), "max_indel_size");
      validateInteger(maxSvBreakpointDistance(), "max_sv_breakpoint_distance");
      validateInteger(maxVariantLengthDifference(), "max_variant_length_difference");
      validateInteger(rescueWindowSize(), "rescue_window_size");
      validateCallsetOrVcf(predictedCallset(), predictedVcf(), "predicted", true);
      validateCallsetOrVcf(trueCallset(), trueVcf(), "true", true);
      validateCallsetOrVcf(knownFalsePositivesCallset(), knownFalsePositivesVcf(), "true", false);
      return this;
    }

    private void validateCallsetOrVcf(
        Optional<String> callset, Optional<File> vcf, String type, boolean required) {
      boolean fileIsPresent = vcf.isPresent();
      if (callset.isPresent()) {
        if (fileIsPresent) {
          throw new IllegalArgumentException(
              String.format("Cannot specifiy both --%s_callset and --%s_vcf", type));
        }
      } else if (!fileIsPresent && required) {
        throw new IllegalArgumentException(
            String.format("Must specifiy either --%s_callset or --%s_vcf", type));
      }
    }

    private void validateDouble(Optional<Double> optional, String flag) {
      if (optional.isPresent()) {
        double d = optional.get();
        if (d < 0.0 || 1.0 < d) {
          throw new IllegalArgumentException(
              String.format("--%s must be between 0.0 and 1.0", flag));
        }
      }
    }

    private void validateInteger(Optional<Integer> optional, String flag) {
      if (optional.isPresent() && optional.get() < 0) {
        throw new IllegalArgumentException(String.format("--%s must be nonnegative", flag));
      }
    }
  }

  public static class FileValidator implements IValueValidator<File> {

    private abstract static class Test {

      private final String reason;

      Test(String reason) {
        this.reason = reason;
      }

      abstract boolean test(File file);

      final void test(String flag, File file) throws ParameterException {
        if (!test(file)) {
          throw new ParameterException(
              String.format("File \"%s\" specified in the \"%s\" flag %s", file, flag, reason));
        }
      }
    }

    private static final Test[] TESTS = {
        new Test("doesn't exist") {
          @Override boolean test(File file) {
            return file.exists();
          }
        },
        new Test("is not a file") {
          @Override boolean test(File file) {
            return file.isFile();
          }
        },
        new Test("is not readable") {
          @Override boolean test(File file) {
            return file.canRead();
          }
        }
    };

    @Override public void validate(String flag, File file) throws ParameterException {
      for (Test test : TESTS) {
        test.test(flag, file);
      }
    }
  }

  @Parameters(separators = "=")
  public static class MainCommand {

    @Parameter(names = { "--lhs" })
    private String leftHandSide;

    @Parameter(names = { "--rhs" })
    private String rightHandSide;

    public Optional<String> leftHandSide() {
      return Optional.fromNullable(leftHandSide);
    }

    public Optional<String> rightHandSide() {
      return Optional.fromNullable(rightHandSide);
    }
  }

  @Parameters(separators = "=")
  public static class NormalizeCommand {

    @Parameter(
        names = { "--callset" })
    private String callset;

    @Parameter(
        names = { "--clean_only" })
    private boolean cleanOnly;

    @Parameter(
        names = { "--max_indel_size" })
    private Integer maxIndelSize;

    @Parameter(
        names = { "--out" },
        converter = FileConverter.class)
    private File out;

    @Parameter(
        names = { "--reference_fasta" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File referenceFasta;

    @Parameter(
        names = { "--reference_fasta_index" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File referenceFastaIndex;

    @Parameter(
        names = { "--vcf" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File vcf;

    public Optional<String> callset() {
      return Optional.fromNullable(callset);
    }

    public boolean cleanOnly() {
      return cleanOnly;
    }

    public Optional<Integer> maxIndexSize() {
      return Optional.fromNullable(maxIndelSize);
    }

    public Optional<File> out() {
      return Optional.fromNullable(out);
    }

    public File referenceFasta() {
      return referenceFasta;
    }

    public Optional<File> referenceFastaIndex() {
      return Optional.fromNullable(referenceFastaIndex);
    }

    NormalizeCommand validate() {
      if (null == referenceFasta) {
        throw new IllegalArgumentException(
            "Flag \"--reference_fasta\" is required");
      }
      if (null == callset && null == vcf) {
        throw new IllegalArgumentException(
            "One of \"--callset\" or \"--vcf\" is required");
      }
      if (null != callset && null != vcf) {
        throw new IllegalArgumentException(
            "Only one of \"--callset\" or \"--vcf\" can be specified");
      }
      if (null == out) {
        throw new IllegalArgumentException(
            "Flag \"--out\" is required");
      }
      return this;
    }

    public Optional<File> vcf() {
      return Optional.fromNullable(vcf);
    }
  }

  @Parameters(separators = "=")
  public static class SetPrefsCommand {

    public enum AuthorizationMethod {
      API_KEY,
      CLIENT_SECRETS,
      SERVICE_ACCOUNT
    }

    public static class AuthorizationMethodConverter
        implements IStringConverter<AuthorizationMethod> {
      @Override public AuthorizationMethod convert(String name) {
        return AuthorizationMethod.valueOf(name);
      }
    }

    @Parameter(
        names = { "--apiKey" })
    private String apiKey;

    @Parameter(
        names = { "--authorizationMethod" },
        converter = AuthorizationMethodConverter.class)
    private AuthorizationMethod authorizationMethod;

    @Parameter(
        names = { "--clientSecretsFile" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File clientSecretsFile;

    @Parameter(
        names = { "--serviceAccountId" })
    private String serviceAccountId;

    @Parameter(
        names = { "--serviceAccountP12File" },
        converter = FileConverter.class,
        validateValueWith = FileValidator.class)
    private File serviceAccountP12File;

    public Optional<String> apiKey() {
      return Optional.fromNullable(apiKey);
    }

    public Optional<AuthorizationMethod> authorizationMethod() {
      return Optional.fromNullable(authorizationMethod);
    }

    public Optional<File> clientSecretsFile() {
      return Optional.fromNullable(clientSecretsFile);
    }

    public Optional<String> serviceAccountId() {
      return Optional.fromNullable(serviceAccountId);
    }

    public Optional<File> serviceAccountP12File() {
      return Optional.fromNullable(serviceAccountP12File);
    }
  }

  public static class ShowPrefsCommand {

    @Parameter(
        names = { "--apiKey" })
    private boolean apiKey;

    @Parameter(
        names = { "--authorizationMethod" })
    private boolean authorizationMethod;

    @Parameter(
        names = { "--clientSecretsFile" })
    private boolean clientSecretsFile;

    @Parameter(
        names = { "--serviceAccountId" })
    private boolean serviceAccountId;

    @Parameter(
        names = { "--serviceAccountP12File" })
    private boolean serviceAccountP12File;

    public boolean apiKey() {
      return apiKey;
    }

    public boolean authorizationMethod() {
      return authorizationMethod;
    }

    public boolean clientSecretsFile() {
      return clientSecretsFile;
    }

    public boolean serviceAccountId() {
      return serviceAccountId;
    }

    public boolean serviceAccountP12File() {
      return serviceAccountP12File;
    }
  }

  protected abstract void diff(DiffCommand diff) throws Exception;
  protected abstract void main(MainCommand command) throws Exception;
  protected abstract void normalize(NormalizeCommand command) throws Exception;
  public final void parse(String... args) throws Exception {
    DiffCommand diff = new DiffCommand();
    MainCommand main = new MainCommand();
    NormalizeCommand normalize = new NormalizeCommand();
    SetPrefsCommand setPrefs = new SetPrefsCommand();
    ShowPrefsCommand showPrefs = new ShowPrefsCommand();
    JCommander jCommander = new JCommander(main);
    jCommander.setProgramName("smash4j");
    jCommander.addCommand("diff", diff);
    jCommander.addCommand("normalize", normalize);
    jCommander.addCommand("setprefs", setPrefs);
    jCommander.addCommand("showprefs", showPrefs);
    try {
      jCommander.parse(args);
      String command = jCommander.getParsedCommand();
      if (null != command) {
        switch (command) {
          case "diff":
            diff(diff.validate());
            break;
          case "normalize":
            normalize(normalize.validate());
            break;
          case "setprefs":
            setPrefs(setPrefs);
            break;
          case "showprefs":
            showPrefs(showPrefs);
            break;
          default:
            throw new IllegalStateException(String.format("Unrecognized command: \"%s\"", command));
        }
      } else {
        main(main);
      }
    } catch (MissingCommandException e) {
      main(main);
    }
  }
  protected abstract void setPrefs(SetPrefsCommand command) throws Exception;

  protected abstract void showPrefs(ShowPrefsCommand command) throws Exception;
}

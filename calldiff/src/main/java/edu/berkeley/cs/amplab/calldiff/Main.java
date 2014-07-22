package edu.berkeley.cs.amplab.calldiff;

import com.google.api.services.genomics.Genomics;

import java.io.File;
import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Main {

  private static final Comparator<Call> COMPARATOR = Comparator
      .comparing(Call::contig)
      .thenComparing(Call::position);

  private static CallScanner callScanner(String name, CommandLine commandLine,
      Optional<String> vcfFile, Optional<String> sampleId, Optional<String> callsetId)
      throws GeneralSecurityException, IOException {
    boolean
        useVcfFile = vcfFile.isPresent(),
        useCallset = callsetId.isPresent();
    if (!useVcfFile && !useCallset) {
      throw new IllegalStateException(
          String.format("Specify one of --%s_vcf or --%s_callset_id", name, name));
    } else if (!useVcfFile && !useCallset) {
      File file = new File(vcfFile.get());
      return sampleId.map(sample -> VcfCallScanner.create(file, sample))
          .orElse(VcfCallScanner.create(file));
    } else if (!useVcfFile && !useCallset) {
      return ApiCallScanner.create(createGenomics(commandLine.apiKey(),
          commandLine.clientSecretsFile(), commandLine.serviceAccountId(), commandLine.p12File(),
          commandLine.rootUrl(), commandLine.timeout()), callsetId.get());
    }
    throw new IllegalStateException(
        String.format("Specify only one of --%s_vcf or --%s_callset_id", name, name));
  }

  private static Genomics createGenomics(
      Optional<String> apiKey,
      Optional<String> clientSecretsFile,
      Optional<String> serviceAccountId,
      Optional<String> p12File,
      Optional<String> rootUrl,
      Optional<Integer> timeout) throws GeneralSecurityException, IOException {
    boolean useApiKey = apiKey.isPresent(),
        useClientSecrets = clientSecretsFile.isPresent(),
        useServiceAccount = serviceAccountId.isPresent() && p12File.isPresent();
    GenomicsFactory.Builder builder = GenomicsFactory.builder();
    rootUrl.ifPresent(url -> builder.setRootUrl(url));
    timeout.ifPresent(ms -> builder.setConnectTimeout(ms).setReadTimeout(ms));
    GenomicsFactory factory = builder.build();
    if (!useApiKey && !useClientSecrets && !useServiceAccount) {
      throw new IllegalStateException("Specify one of { --api_key, --client_secrets_file, "
          + "(--service_account_id, --p12_file) }");
    } else if (useApiKey && !useClientSecrets && !useServiceAccount) {
      return factory.fromApiKey(apiKey.get());
    } else if (!useApiKey && useClientSecrets && !useServiceAccount) {
      return factory.fromClientSecretsFile(new File(clientSecretsFile.get()));
    } else if (!useApiKey && !useClientSecrets && useServiceAccount) {
      return factory.fromServiceAccount(serviceAccountId.get(), new File(p12File.get()));
    }
    throw new IllegalStateException("Specify only one of { --api_key, --client_secrets_file, "
        + "(--service_account_id, --p12_file) }");
  }

  private static FastaReader fastaReader(Optional<String> referenceFasta,
      Optional<String> referenceFai) throws Exception {
    try {
      return referenceFasta.map(File::new)
          .map(fastaFile -> {
                try {
                  return referenceFai.map(File::new)
                      .map(faiFile -> {
                            try {
                              return FastaReader.create(fastaFile, faiFile);
                            } catch (IOException e) {
                              throw ExceptionWrapper.wrap(e);
                            }
                          })
                      .orElse(FastaReader.create(fastaFile));
                } catch (IOException e) {
                  throw ExceptionWrapper.wrap(e);
                }
              })
          .orElseThrow(() -> new IllegalArgumentException("--reference_fasta is required"));
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }

  public static void main(String[] args) throws Exception {
    try {
      CommandLine.parse(args).ifPresent(commandLine -> {
        try {
          System.out.println(fastaReader(commandLine.referenceFasta(), commandLine.referenceFai())
              .read(reference -> {
                try {
                  return callScanner(
                          "lhs",
                          commandLine,
                          commandLine.lhsVcf(),
                          commandLine.lhsSampleId(),
                          commandLine.lhsCallsetId())
                      .scan(lhs -> {
                        try {
                          boolean presorted = commandLine.presorted();
                          return callScanner(
                                  "rhs",
                                  commandLine,
                                  commandLine.rhsVcf(),
                                  commandLine.rhsSampleId(),
                                  commandLine.rhsCallsetId())
                              .scan(rhs -> OutputTuple
                                  .calldiff(
                                      reference,
                                      presorted ? lhs : sort(lhs),
                                      presorted ? rhs : sort(rhs))
                                  .collect(DiffStats.builder()));
                        } catch (GeneralSecurityException | IOException e) {
                          throw ExceptionWrapper.wrap(e);
                        }
                      });
                } catch (GeneralSecurityException | IOException e) {
                  throw ExceptionWrapper.wrap(e);
                }
              }));
        } catch (Exception e) {
          throw ExceptionWrapper.wrap(e);
        }
      });
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }

  private static Stream<Call> sort(Stream<Call> stream) {
    return sort(stream, COMPARATOR);
  }

  private static <X> Stream<X> sort(Stream<X> stream, Comparator<? super X> comparator) {
    List<X> list = stream.collect(Collectors.toCollection(ArrayList::new));
    Collections.sort(list, comparator);
    return list.stream();
  }
}

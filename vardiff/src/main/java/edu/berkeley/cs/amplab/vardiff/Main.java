package edu.berkeley.cs.amplab.vardiff;

import java.io.File;
import java.io.IOException;
import java.util.Optional;

public class Main {

  private static class ExceptionWrapper extends RuntimeException {

    ExceptionWrapper(IOException cause) {
      super(cause);
    }

    @Override public IOException getCause() {
      return (IOException) super.getCause();
    }
  }

  private static CallScanner
      callScanner(String name, Optional<String> vcfFile, Optional<String> sampleId) {
    return vcfFile.map(File::new)
        .map(
            vcf -> sampleId.map(id -> VcfCallScanner.create(vcf, id))
                .orElse(VcfCallScanner.create(vcf)))
        .orElseThrow(
            () -> new IllegalArgumentException(String.format("--%s_vcf is required", name)));
  }

  private static FastaReader fastaReader(Optional<String> referenceFasta,
      Optional<String> referenceFai) throws IOException {
    try {
      return referenceFasta.map(File::new)
          .map(fastaFile -> {
                try {
                  return referenceFai.map(File::new)
                      .map(faiFile -> {
                            try {
                              return FastaReader.create(fastaFile, faiFile);
                            } catch (IOException e) {
                              throw new ExceptionWrapper(e);
                            }
                          })
                      .orElse(FastaReader.create(fastaFile));
                } catch (IOException e) {
                  throw new ExceptionWrapper(e);
                }
              })
          .orElseThrow(() -> new IllegalArgumentException("--reference_fasta is required"));
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }

  public static void main(String[] args) throws IOException {
    try {
      CommandLine commandLine = CommandLine.parse(args);
      System.out.println(
          fastaReader(commandLine.referenceFasta(), commandLine.referenceFai())
              .read(reference -> {
                try {
                  return callScanner("lhs", commandLine.lhsVcf(), commandLine.lhsSampleId())
                      .scan(lhs -> {
                        try {
                          return callScanner("rhs", commandLine.rhsVcf(), commandLine.rhsSampleId())
                              .scan(rhs -> OutputTuple.vardiff(reference, lhs, rhs)
                                  .collect(DiffStats.builder()));
                        } catch (IOException e) {
                          throw new ExceptionWrapper(e);
                        }
                      });
                } catch (IOException e) {
                  throw new ExceptionWrapper(e);
                }
              }));
    } catch (ExceptionWrapper e) {
      throw e.getCause();
    }
  }
}

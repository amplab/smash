package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import edu.berkeley.cs.amplab.vcfparser.Header;
import edu.berkeley.cs.amplab.vcfparser.MetaInformation;
import edu.berkeley.cs.amplab.vcfparser.VcfReader;
import edu.berkeley.cs.amplab.vcfparser.VcfRecord;
import edu.berkeley.cs.amplab.vcfparser.VcfWriter;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

public class Main {

  /**
   * Right now, all this does is copy the VCF data in the file specified by the
   * --in command-line flag to the file specified by the --out command-line flag.
   */
  public static void main(String[] args) throws Exception {
    Optional<CommandLine> optional = CommandLine.parse(args);
    if (optional.isPresent()) {
      final CommandLine commandLine = optional.get();
      try (InputStream in = new FileInputStream(commandLine.in())) {
        try (OutputStream out = new FileOutputStream(commandLine.out())) {
          VcfReader.from(in).read(
              new VcfReader.Callback<OutputStream>() {
                @Override
                public OutputStream readVcf(
                    MetaInformation metaInformation,
                    Header header,
                    FluentIterable<VcfRecord> records) {
                  return VcfWriter.to(out).write(metaInformation, header, records);
                }
              });
        }
      }
    }
  }
}

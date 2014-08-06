NOTE
================
The Java version of SMaSH is not yet officially supported and should be considered experimental.


Before you build
================

You need to have [JDK 1.8](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
and [Apache Maven](http://maven.apache.org/download.cgi) installed in order to build the code.

How to build
============

To build the code, just run `mvn assembly:assembly` from within the `calldiff` directory.
This should produce the fat jar file `target/calldiff-jar-with-dependencies.jar`.

Command line usage
===================

Run `java -jar target/calldiff-jar-with-dependencies.jar --help`. You should see output
similar to the following:

    Usage: java -jar calldiff-jar-with-dependencies.jar [options]
      Options:
      --api_key             The API key used to authenticate to your Google Cloud
                            project
      --client_secrets_file The client secrets file used to authorize access to your
                            Google Cloud project
      --help                Print the help message
      --lhs_callset_id      The callset id to use on the left hand side of the
                            comparison
      --lhs_sample_id       The sample id to use on the left hand side of the
                            comparison
      --lhs_vcf             The path to the VCF file to use on the left hand side of
                            the comparison
      --p12_file            The P12 file containing the private key that authorizes
                            the service account for your Google Cloud Project
      --presorted           Skip sorting the input because it is already properly
                            sorted
      --reference_fai       The FASTA index file for the reference sequence
      --reference_fasta     The FASTA file for the reference sequence
      --rhs_callset_id      The callset id to use on the right hand side of the
                            comparison
      --rhs_sample_id       The sample id to use on the right hand side of the
                            comparison
      --rhs_vcf             The path to the VCF file to use on the right hand side
                            of the comparison
      --root_url            The URL to communicate with to fetch variants from the
                            cloud
      --service_account_id  The email address for the service account used to
                            authorize your Google Cloud project
      --timeout             The connect and read timeouts to use when making
                            requests to the cloud


The options starting with `--lhs` refer to the callset on the left-hand side of the
comparison, and the options starting with `--rhs` refer to the callset on the right
hand side of the comparison. Each comparison takes exactly one callset on the left
and one on the right. Callsets may originate either as columns of a VCF file or as
data served from the cloud. To specify that the callset comes from the cloud, using
the `--[lr]_callset_id` flag. To specify a column from a VCF file, use the
`--[lr]hs_vcf` flag. If the VCF file has more than one sample in it, you must
specify the sample to use with the `--[lr]hs_sample_id` flag.

All comparisons require a reference sequence, supplied via a FASTA file using the
`--reference_fasta` flag. The code will attempt to find a FASTA index file for the
given reference by looking for a file with the exact same name as the reference
file, except with a `.fai` suffix. If no FASTA index file for the reference is
found, the index will be precomputed in advance of making the comparison. You can
always explicitly tell the program where the index file is using the
`--reference_fai` flag.

If you are comparing callsets that are served from the cloud, you must also provide
command line flags for specifying the authentication mechanism. Right now, there
are 3 supported ways to authenticate: Using an API Key (the `--api_key` flag),
using client secrets (The `--client_secrets_file` flag) or using a service
account (The `--service_account_id` and `--p12_file` flags).


The code
========

See [this doc](https://docs.google.com/document/d/1gEpZVsNgkZAjbgudI-KsRILXi7TSD459a4LMAP9s1LQ/edit?usp=sharing)
for more details on the methodology behind this code.

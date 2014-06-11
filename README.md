SMaSH
=====

A benchmarking toolkit for variant calling. For a full description of the SMaSH toolkit, including our datasets, results, and instructions for running SMaSH on EC2, see our website at smash.cs.berkeley.edu.

If you have any questions, please visit our support forum at https://groups.google.com/forum/#!forum/smash-benchmarking.

SMaSH Overview
--------------
The simplest SMaSH use case is to compare the VCF output of a variant calling algorithm against a validated truth set.
We do this by calling the bench.py script with three main arguments: the ground truth VCF, the predicted VCF, and a reference FASTA.

## Evaluation

SMaSH evaluates three main types of variants: SNPs, indels, and structural variants. Indels are defined as a variant shorter than 50
base pairs; longer variants are considered structural variants. Indels are further subdivided into insertions, deletions, inversions, and other; structural variants are likewise broken into insertions, deletions, and others.

First, SMaSH checks all locations where both the true VCF and predicted VCF have variants. SNPs and indels are considered matches if the reference allele and first alternative allele are the same. Structural variants are considered matches if they are the same type of variant (insertion, deletion, or many-to-many) and if their position is with a certain breakpoint tolerance of the ground truth variant. (SV tolerance can be specified when calling SMaSH; the default value is 50 base pairs).

However, the same underlying variants can be represented in different ways in the VCF format. SMaSH addresses this ambiguity through its rescue algorithm: for each variant in the truth VCF not matched in the predicted VCF, SMaSH expands the sequence described by the true and predicted callsets to a defined window size (can be specified at runtime, defaults to 100bp). If the sequences match, the true callset variants are marked as true positives, and the predicted callset variants are removed from the set of false positives.

## Outputs

After evaluation, SMaSH reports results as human-readable text (default) or as tab-delimited output. The results include 
counts of true positives, false positives, and false negatives; as well as precision and recall calculations with bounded error based on the error rates passed in at runtime; and some genotyping accuracy metrics. SMaSH can also output a fully annotated VCF file containing all variants from both the true and predicted VCF input files, tagged as true positive, false positive, false negative, or rescued.

## Known False Positive Mode

For cases in which we do not have validated truth calls for the entire genome, SMaSH may be called with an additional known true positive VCF. With this option, SMaSH will further mark predicted calls as conflicting with known false positives if they are at the same location, of the same variant type, and have the same ref sequence. Precision will also be calculated with this false positive count.

## Normalization

SMaSH also includes a normalization script which first cleans the VCF by upper-casing all alleles, removing homozygous reference and monomorphic allelic calls, and so on. The script then left-normalizes all calls as much as possible as a step toward standardizing ambiguous calls in a repetitive sequence. Although left-normalization is not strictly necessary to run SMaSH, it improves accuracy and reduces runtime by reducing the number of times the rescue algorithm needs to be called. Normalization can be performed as part of the evaluation script call by adding the --normalize option.

## Assumptions and Limitations

* At present SMaSH does not evaluate phasing.
* SMaSH does not handle compound heterozygous variants.


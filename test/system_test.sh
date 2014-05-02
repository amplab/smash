#!/usr/bin/env sh

# Directory definitions
test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
smash_dir=$test_dir/..
smashbenchmarking_dir=$smash_dir/smashbenchmarking
scripts_dir=$smash_dir/scripts
tmp_dir=`mktemp -d`

# Normalize inputs
for vcf_file in $test_dir/system_test_true.vcf $test_dir/system_test_pred.vcf
do
  head -n 1000 $vcf_file | grep ^# > $tmp_dir/hdr
  cat $vcf_file | python $smashbenchmarking_dir/normalize_vcf.py \
      $test_dir/ref.fasta myvcf 50 2> /dev/null | grep -v ^# | sort -nk2,2 | \
      perl $scripts_dir/sortByRef.pl - $test_dir/ref.fasta.fai | cat \
      $tmp_dir/hdr - > $tmp_dir/$(basename $vcf_file).normalized
done

# Run the benchmarker
/usr/bin/env python $smashbenchmarking_dir/bench.py \
    $tmp_dir/system_test_true.vcf.normalized \
    $tmp_dir/system_test_pred.vcf.normalized $test_dir/ref.fasta --snp_err 0.0 \
    --indel_err 0.0 --sv_err 0.0 --sv_bp 100 -w 50 1> $tmp_dir/output \
    2> /dev/null

# Diff against golden output
diff $test_dir/system_test.golden_output $tmp_dir/output > $tmp_dir/diff_output
exit_code=$?

if [ $exit_code -eq 0 ]; then
  # clean up temp files if the test passed
  rm -rf $tmp_dir
else
  echo SMaSH output differs from golden output. See $tmp_dir/diff_output to \
      see why.
fi

# Report the 'diff' exit code as this test's exit code
exit $exit_code


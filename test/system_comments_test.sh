#!/usr/bin/env bash

# Directory definitions
test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
smash_dir=$test_dir/..
smashbenchmarking_dir=$smash_dir/smashbenchmarking
scripts_dir=$smash_dir/scripts
tmp_dir=`mktemp -d /tmp/systest.XXXX`

# strip comments out of test vcf files
grep -v '^--' $test_dir/true_comments.vcf > $tmp_dir/true_comments.vcf
grep -v '^--' $test_dir/pred_comments.vcf > $tmp_dir/pred_comments.vcf

# Run the benchmarker
/usr/bin/env python $smashbenchmarking_dir/bench.py \
    $tmp_dir/true_comments.vcf \
    $tmp_dir/pred_comments.vcf $test_dir/ref.fasta $test_dir/ref.fasta.fai --snp_err 0.0 \
    --indel_err 0.0 --sv_err 0.0 --sv_bp 100 -w 50 --output tsv --normalize --output_vcf $tmp_dir/output_vcf 1> $tmp_dir/output \
    2> /dev/null

# Diff against golden output
diff $test_dir/system_comments_test.golden_output $tmp_dir/output --ignore-space-change --ignore-matching-lines=^# > $tmp_dir/diff_output
exit_code=$?

diff $test_dir/system_comments_golden_output.vcf $tmp_dir/output_vcf --ignore-space-change --ignore-matching-lines=^# > $tmp_dir/diff_output_vcf 
vcf_exit_code=$?

if [ $exit_code -eq 0 ]; then
  if [ $vcf_exit_code -eq 0 ]; then
    # both parts of the test passed!
    # clean up temp files
    rm -rf $tmp_dir
  else
    echo SMaSH vcf output differs from golden_output. See $tmp_dir/diff_output_vcf to see why.
  fi
else
  echo SMaSH output differs from golden output. See $tmp_dir/diff_output to \
      see why.
fi

# Report the 'diff' exit code as this test's exit code
exit $exit_code


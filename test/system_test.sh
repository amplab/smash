#!/usr/bin/env bash

# Directory definitions
test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
smash_dir=$test_dir/..
smashbenchmarking_dir=$smash_dir/smashbenchmarking
scripts_dir=$smash_dir/scripts
tmp_dir=`mktemp -d /tmp/systest.XXXX`

# Run the benchmarker
/usr/bin/env python $smashbenchmarking_dir/bench.py \
    $test_dir/system_test_true.vcf \
    $test_dir/system_test_pred.vcf $test_dir/ref.fasta --snp_err 0.0 \
    --indel_err 0.0 --sv_err 0.0 --sv_bp 100 -w 50 --normalize 1> $tmp_dir/output \
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


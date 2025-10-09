#!/bin/bash


set -e  # Exit on error

echo "Running seq_merge..."
./seq_merge --in1 test/testing_reads.fastq --bc-start 0 --bc-len 15 --out1 test/test_output.fastq --threads 10

echo "Comparing output to target..."
result=$(python3 test/compare_fastqs.py test/test_output.fastq test/target_output.fastq)

if [ "$result" = "True" ]; then
    echo "SUCCESS: Output matches target. seq_merge is functioning correctly."
    rm test/test_output.fastq
    exit 0
else
    echo "FAILURE: Output does not match target. seq_merge output differs from expected."
    exit 1
fi

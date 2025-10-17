#!/bin/bash


set -e  # Exit on error

echo "Running bcar_sort..."
./bcar_sort --in test/testing_reads.fastq --bc-start 0 --bc-len 15 --out test/test_sorted.fastq --temp .
echo "Running bcar_merge..."
./bcar_merge --in test/test_sorted.fastq --bc-start 0 --bc-len 15 --out test/test_output.fastq --threads 10

echo "Comparing output to target..."
result=$(python3 test/compare_fastqs.py test/test_output.fastq test/target_output.fastq)

if [ "$result" = "True" ]; then
    echo "SUCCESS: Output matches target. bcar is functioning as expected."
    rm test/test_output.fastq
    rm test/test_sorted.fastq
    exit 0
else
    echo "FAILURE: Output does not match target. bcar output differs from expected."
    exit 1
fi

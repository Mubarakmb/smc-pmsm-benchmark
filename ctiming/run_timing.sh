#!/usr/bin/env bash
# Compile and time the six control laws (3 repetitions, 2e7 calls each; FOSMC 4e5 calls).
# Results are machine-dependent; the paper reports the median of 3 runs on an Intel Xeon @ 2.8 GHz.
set -euo pipefail
cd "$(dirname "$0")"
gcc -O2 -std=c11 ctrl_bench.c -lm -o ctrl_bench
for r in 1 2 3; do
  echo "# run $r"
  ./ctrl_bench 20000000
done | tee timing_output.txt

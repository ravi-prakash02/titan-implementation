#!/bin/bash

echo "benchmarking Titan ..."
for t in 4 8 12 16; do
    tmpfile="titan_${t}t_raw.log"

    taskset -c 0-$((t-1)) env \
    RAYON_NUM_THREADS=$t \
    RUST_TEST_THREADS=1 \
    cargo test titan_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
    > "$tmpfile"
    grep '^[0-9]' "$tmpfile" > titan_${t}t_clean.csv
    rm "$tmpfile"
done
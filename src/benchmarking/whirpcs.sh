#!/bin/bash

# COPY THIS FILE INSIDE THE MODIFIED WHIR REPO AND RUN IT THERE (MODIFIED: prints commit and prove time separately)

echo "m,commit_ms,prove_ms,verifier_ms,proof_bytes" > whir.csv

for m in 18 20 22 24 26; do
  echo "Running WHIR m=$m (5 repetitions)"

  commit_sum=0
  prove_sum=0
  verifier_sum=0
  size_sum=0
  reps=5

  for ((i=1; i<=reps; i++)); do
    echo "  Run $i/$reps"

    out=$(taskset -c 0 env RAYON_NUM_THREADS=1 cargo run --release -- \
        --type PCS \
        --sec ProvableList \
        --num-variables $m \
        --pow-bits 0 \
        --field Field256 \
        --evaluations 1 \
        --rate 1
    )

  
    # Convert time string to ms
    convert_to_ms() {
      val=$(echo "$1" | grep -oE '[0-9.]+')
      unit=$(echo "$1" | grep -oE '[a-zA-Zµ]+')
      case "$unit" in
        s)  echo "$(awk "BEGIN {print $val*1000}")" ;;
        ms) echo "$val" ;;
        µs) echo "$(awk "BEGIN {print $val/1000}")" ;;
        *)  echo "$val" ;;
      esac
    }

    # Extract raw values
    commit_raw=$(echo "$out" | grep "Commit time" | awk '{print $3}')
    prove_raw=$(echo "$out" | grep "Prove time" | awk '{print $3}')
    verifier_raw=$(echo "$out" | grep "Verifier time" | awk '{print $3}')
    size_raw=$(echo "$out" | grep "Proof size" | awk '{print $3}')

    # Convert to numeric values
    commit_ms=$(convert_to_ms "$commit_raw")
    prove_ms=$(convert_to_ms "$prove_raw")
    verifier_ms=$(convert_to_ms "$verifier_raw")
    proof_bytes=$(echo "$size_raw" | awk '{print int($1*1024)}')  # KiB to bytes

    commit_sum=$(awk "BEGIN {print $commit_sum+$commit_ms}")
    prove_sum=$(awk "BEGIN {print $prove_sum+$prove_ms}")
    verifier_sum=$(awk "BEGIN {print $verifier_sum+$verifier_ms}")
    size_sum=$(awk "BEGIN {print $size_sum+$proof_bytes}")
  done

  # Compute averages
  commit_avg=$(awk "BEGIN {print $commit_sum/$reps}")
  prove_avg=$(awk "BEGIN {print $prove_sum/$reps}")
  verifier_avg=$(awk "BEGIN {print $verifier_sum/$reps}")
  size_avg=$(awk "BEGIN {print $size_sum/$reps}")

  echo "$m,$commit_avg,$prove_avg,$verifier_avg,$size_avg" >> whir.csv
done

# Titan: Commitment Scheme Benchmarking (Rust)

This project benchmarks **Titan**, our polynomial commitment scheme, against existing schemes.

The evaluation compares performance across the following metrics:

- Commit time  
- Evaluation (proof generation) time  
- Verification time  
- Proof size  

---

## Benchmarking Setup

All experiments are executed under strict conditions to ensure fair comparison:

- Single-core execution using `taskset` (Linux only)
- Parallelism disabled:
  - `RAYON_NUM_THREADS=1`
  - `RUST_TEST_THREADS=1`
- Tests run in release mode
- CPU affinity is fixed where supported

---

## Platform Compatibility

The benchmarking commands use `taskset`, which is **only available on Linux systems** (via `util-linux`).

---

## Running Experiments

### Titan (Proposed Scheme)

```bash
taskset -c 0 env RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 \
cargo test titan_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
| grep '^[0-9]' > titan_clean.csv
```

Generate plots:

```bash
python titan_plot.py
```

---

### Dory (Baseline)

```bash
taskset -c 0 env RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 \
cargo test dory_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
| grep '^[0-9]' > dory_clean.csv
```

---

## Output Format

The cleaned CSV files (`*_clean.csv`) contain numeric benchmark results. Each row follows:

```text
input_size, commit_time, eval_time, verify_time, proof_size
```

---


## Security Disclaimer

**Titan is an experimental commitment scheme and has not undergone formal security analysis or auditing.**

- It may contain **undiscovered vulnerabilities**
- It should **NOT be used in production systems**
- It is intended **strictly for academic and research purposes**

---
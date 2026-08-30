# Titan Benchmarking

This project implements Titan, our polynomial commitment scheme (PCS) and benchmarks it against some of the existing schemes. 

The PCS evaluation compares performance across the following metrics:

- Commit time  
- Evaluation (proof generation) time  
- Verification time  
- Proof size  

This project also implements the MicroSpartan PIOP with Titan PCS to obtain a transparent SNARK.

---

## Benchmarking Setup

All experiments are executed under strict conditions to ensure fair comparison:

- Single-core execution using `taskset` (Linux only, see [Platform Compatibility](#platform-compatibility))
- Parallelism disabled using:
  - `RAYON_NUM_THREADS=1`
  - `RUST_TEST_THREADS=1`
- Tests run in release mode
- CPU affinity is fixed where supported

---
## Installing Rust

This project is implemented in Rust. To build and run the benchmarks, you need to install the Rust toolchain. For detailed installation instructions, refer to the [official Rust documentation](https://doc.rust-lang.org/book/ch01-01-installation.html).

### Install via rustup (recommended)

Run the following command in a terminal:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then follow the on-screen instructions. After installation, restart your shell (or close and reopen your terminal).

---

### Verify Installation

Check that Rust is installed correctly:

```bash
rustc --version
cargo --version
```

---

### Updating Rust

To update to the latest version:

```bash
rustup update
```

---

### Platform Notes for Rust

- Linux and macOS are fully supported  
- On Windows, use PowerShell or install via WSL for best compatibility  
- For benchmarking Titan PCS, Linux is recommended (see [Platform Compatibility](#platform-compatibility))

---
## Getting Started

Download the repository as a ZIP file, extract it locally, and navigate into the extracted directory using a terminal:

```bash
cd path/to/titan-implementation-BB26
```

All the provided benchmarking and other commands should be run from the root of this directory.

---

## Running Experiments

To benchmark Titan PCS, run:

```bash
taskset -c 0 env RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 \
cargo test titan_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
| grep '^[0-9]' > src/benchmarking/titan_clean.csv
```

To generate plots, Python 3 must be installed. Check whether Python is installed:

```bash
python3 --version
```

If Python 3 is not installed, install it using your system's package manager. For example, on Ubuntu/Debian:

```bash
sudo apt update
sudo apt install python3
```

Then generate the plots:

```bash
python3 src/benchmarking/titan_plot.py
```

---
### Other Schemes (Baselines)

Benchmarking results for other schemes (Brakedown, Dory, Hyrax, Kopis, Whir) have already been included. Follow the steps below to reproduce them.

To reproduce the results of Dory and Kopis, run their corresponding scaling experiments (using the **same command structure and environment settings** as above). 

Example, for Dory: 
 ```bash
taskset -c 0 env RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 \
cargo test dory_BLS_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
| grep '^[0-9]' > src/benchmarking/dory_BLS_clean.csv
```

- **Brakedown & Hyrax:** Benchmarks were obtained by modifying the [`poly-commit`](https://github.com/arkworks-rs/poly-commit) repository. Follow the instructions in `src/benchmarking/brakedown_hyrax_benchmarking_steps.txt` file to reproduce the results.

- **Whir:**
Benchmarks were obtained by modifying the [`Whir`](https://github.com/WizardOfMenlo/whir.git) repository. To reproduce them:
  - Make sure [Git](https://git-scm.com/) is installed:
     ```bash
     git --version
     ```
     If Git is not installed, follow the [Git installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

  - Clone the Whir repository:
     ```bash
     git clone https://github.com/WizardOfMenlo/whir.git
     ```

  - Copy the `whirpcs.sh` benchmarking script (which is present inside `src/benchmarking`) into the Whir repository:
     ```bash
     cp src/benchmarking/whirpcs.sh whir/
     ```

  - Navigate to the Whir repository:
     ```bash
     cd whir
     ```

  - Run the benchmarking script:
     ```bash
     ./whirpcs.sh
     ```
---

## Output Format

The cleaned CSV files (`*_clean.csv`) and the PNG files will be saved in `src/benchmarking`. Each row of the CSV file contains:

```text
input_size, commit_time, eval_time, verify_time, proof_size
```

---
## Platform Compatibility

The benchmarking commands use `taskset`, which is **only available on Linux systems** (via `util-linux`).

On unsupported platforms, you can still run the benchmarks by **removing `taskset -c 0`**, but results may be less consistent due to OS scheduling across multiple CPU cores.

Example (macOS / Windows alternative):

```bash
RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 \
cargo test titan_scaling_experiment --release -- --ignored --nocapture --test-threads=1 \
| grep '^[0-9]' > src/benchmarking/titan_clean.csv
```
---
## Executing MicroSpartan PIOP with Titan PCS

Run:

```bash
cargo test test_spartan_titan_index --release -- --nocapture
```

This test prints:

- Setup time  
- Indexing time  
- Proving time  
- Verification time  
- Proof size  

### Notes

- This test is intended for **functional validation** only and not for strict benchmarking  
- Parameters (e.g., constraint size, number of queries) can be modified directly in the test code (in `src/protocols/spartan_titan.rs`)  
- No special environment configuration is required  

---

## Security Disclaimer

**Titan and MicroSpartan PIOP with Titan PCS are experimental constructions and have not undergone formal security analysis or auditing.**

- It may contain **undiscovered vulnerabilities**
- It should **NOT be used in production systems**
- It is intended **strictly for academic and research purposes**

---

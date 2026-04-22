import pandas as pd
import matplotlib.pyplot as plt

threads = [4, 8, 12, 16]
files = {t: f"titan_{t}t_clean.csv" for t in threads}

data = {}
for t, f in files.items():
    df = pd.read_csv(
        f,
        header=None,
        names=["m", "commit_ms", "prove_ms", "verify_ms", "proof_bytes"]
    )
    data[t] = df


def plot_metric(metric, ylabel, title, filename):
    plt.figure(figsize=(6, 4))

    for t, df in data.items():
        plt.plot(
            df["m"],
            df[metric],
            marker="o",
            label=f"{t} threads",
        )

    plt.xlabel("Number of variables (m)")
    plt.ylabel(ylabel)
    plt.yscale("log")   # important for scaling plots
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# Commit time
plot_metric(
    metric="commit_ms",
    ylabel="Commit time (ms)",
    title="Titan Commit Time vs Number of Variables",
    filename="titan_commit_time_threads.png",
)

# Prove time
plot_metric(
    metric="prove_ms",
    ylabel="Eval time (ms)",
    title="Titan Prove Time vs Number of Variables",
    filename="titan_eval_time_threads.png",
)

# Verify time
plot_metric(
    metric="verify_ms",
    ylabel="Verify time (ms)",
    title="Titan Verify Time vs Number of Variables",
    filename="titan_verify_time_threads.png",
)

# Proof size
plot_metric(
    metric="proof_bytes",
    ylabel="Proof size (bytes)",
    title="Titan Proof Size vs Number of Variables",
    filename="titan_proof_size_threads.png",
)

